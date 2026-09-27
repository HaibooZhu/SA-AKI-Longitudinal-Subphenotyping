"""Read-only model audit: compare saved AUMC inputs with all calcium features missing.

No training, parameter selection, source-input changes, or manuscript replacement.
New predictions remain private on the server. Aggregate comparison is diagnostic.
"""
import os
os.environ.setdefault('OMP_NUM_THREADS', '4')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
from pathlib import Path
import argparse, hashlib, json, sys, platform
import importlib.metadata as md
import numpy as np
import pandas as pd
from autogluon.tabular import TabularPredictor

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def main(a):
    os.umask(0o077)
    a.out.mkdir(exist_ok=False)
    sys.path.insert(0,str(a.legacy))
    from train import metrics
    rows=[]; checks={}; model_hashes={}
    plan={'scope':'AUMC calcium-missing diagnostic comparison; no retraining or model selection',
          'hours':[6,12,18,24,30,36], 'rule':'Set every calcium_ derived input to NaN, leave all other inputs and labels unchanged',
          'source':str(a.source),'script_sha256':sha(Path(__file__)),
          'status':'Independent audit output; not adopted into manuscript',
          'python':platform.python_version(),'autogluon':md.version('autogluon.tabular')}
    (a.out/'plan.json').write_text(json.dumps(plan,indent=2))
    for hour in plan['hours']:
        root=a.source/'windows'/f'window_{hour:02d}h'
        files=[p for p in (root/'predictor').rglob('*') if p.is_file()]
        files += [root/'data/external.pkl',root/'private_predictions_external.pkl',root/'data/train.pkl',root/'data/internal.pkl']
        hashes={str(p):sha(p) for p in files}
        model_hashes.update(hashes)
        df=pd.read_pickle(root/'data/external.pkl')
        saved=pd.read_pickle(root/'private_predictions_external.pkl')
        pd.testing.assert_frame_equal(df[['stay_id','groupHPD']],saved[['stay_id','groupHPD']])
        assert len(df)==2183 and set(df.dataset)=={'aumcdb'}
        X=df.drop(columns=['stay_id','dataset','groupHPD'])
        y=df.groupHPD.to_numpy(dtype=int)
        ca=[c for c in X if c.startswith('calcium_')]
        expected=3+(hour//6)*(hour//6-1)//2
        assert len(ca)==expected
        changed=X.copy(); changed.loc[:,ca]=np.nan
        pd.testing.assert_frame_equal(X.drop(columns=ca),changed.drop(columns=ca))
        assert changed[ca].isna().all().all()
        predictor=TabularPredictor.load(str(root/'predictor'),verbosity=0)
        model=predictor.model_best
        old=predictor.predict_proba(X,model=model).loc[:,[1,2,3]].to_numpy()
        savedp=saved[['p_DR','p_RR','p_PW']].to_numpy()
        error=float(np.max(np.abs(old-savedp)))
        assert error<1e-7,(hour,error)
        new=predictor.predict_proba(changed,model=model).loc[:,[1,2,3]].to_numpy()
        assert np.isfinite(new).all() and np.allclose(new.sum(axis=1),1,atol=1e-6)
        originally_missing=X[ca].isna().all(axis=1).to_numpy()
        assert np.allclose(new[originally_missing],old[originally_missing],atol=1e-7)
        for name,p in [('saved_input',old),('calcium_missing',new)]:
            rows.append({'hours':hour,'version':name,'n':len(y),**metrics(y,p)})
        private=saved.copy()
        private[['p_DR','p_RR','p_PW']]=new
        private.to_pickle(a.out/f'private_predictions_external_{hour:02d}h_calcium_missing.pkl')
        checks[str(hour)]={'model':model,'features':ca,'feature_count':len(ca),
            'saved_prediction_max_abs_error':error,
            'patients_with_any_calcium':int((~originally_missing).sum()),
            'class_counts':df.groupHPD.value_counts().sort_index().to_dict(),
            'non_calcium_inputs_unchanged':True,
            'mean_abs_probability_change':float(np.mean(np.abs(old-new))),
            'max_abs_probability_change':float(np.max(np.abs(old-new))),
            'argmax_changed':int((old.argmax(1)!=new.argmax(1)).sum())}
        if hour==24:
            internal=pd.read_pickle(root/'data/internal.pkl')
            ips=pd.read_pickle(root/'private_predictions_internal.pkl')
            ip=predictor.predict_proba(internal.drop(columns=['stay_id','dataset','groupHPD']),model=model).loc[:,[1,2,3]].to_numpy()
            assert np.max(np.abs(ip-ips[['p_DR','p_RR','p_PW']].to_numpy()))<1e-7
            checks['24']['internal_predictions_reproduced']=True
            train=pd.read_pickle(root/'data/train.pkl')
            checks['24']['calcium_mean_summary']={co:{'n_nonmissing':int(g.calcium_mean.notna().sum()),'median':float(g.calcium_mean.median())} for co,g in pd.concat([train,df]).groupby('dataset')}
        assert all(sha(Path(p))==h for p,h in hashes.items()),'Source files changed'
        print(json.dumps({'hours':hour,'reproduced':error,'old':rows[-2]['auc_ovo_macro'],'new':rows[-1]['auc_ovo_macro'],'old_DR_PW':rows[-2]['auc_DR_PW'],'new_DR_PW':rows[-1]['auc_DR_PW']},ensure_ascii=False),flush=True)
    pd.DataFrame(rows).to_csv(a.out/'comparison_metrics.csv',index=False)
    (a.out/'verification.json').write_text(json.dumps(checks,indent=2))
    (a.out/'source_sha256.json').write_text(json.dumps(model_hashes,indent=2))
    (a.out/'COMPLETE.json').write_text(json.dumps({'source_files_unchanged':True,'six_windows_reproduced':True,'no_fit_called':True,'manuscript_unchanged':True},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--source',type=Path,default=Path('~/jtim_refit_units_20260922').expanduser())
    p.add_argument('--legacy',type=Path,default=Path('~/jtim_autogluon_20260921').expanduser())
    p.add_argument('--out',type=Path,required=True)
    main(p.parse_args())
