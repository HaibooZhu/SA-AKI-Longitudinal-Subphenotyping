"""Recompute all metrics, paired window changes and stratified bootstrap CIs."""
import argparse,json,time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from train import metrics

HOURS=[6,12,18,24,30,36]
NAMES=['auc_ovo_macro','auc_DR_RR','auc_DR_PW','auc_RR_PW']

def auc(y,s):
    n1=y.sum();n0=len(y)-n1
    return (rankdata(s)[y].sum()-n1*(n1+1)/2)/(n1*n0)

def four_auc(y,p):
    macro=[];pairs=[]
    for a,b in [(1,2),(1,3),(2,3)]:
        keep=np.isin(y,[a,b]);pos=y[keep]==a;pa=p[keep,a-1];pb=p[keep,b-1]
        macro.append((auc(pos,pa)+auc(~pos,pb))/2)
        pairs.append(auc(pos,pa/np.maximum(pa+pb,1e-15)))
    return np.array([np.mean(macro),*pairs])

def analyze(root,draws=1000):
    assert (root/'TRAINING_COMPLETE.json').exists()
    tab=pd.read_csv(root/'metrics.csv');assert len(tab)==12
    cis=[];changes=[];verify={'recomputed_evaluations':[],'bootstrap_draws':draws,
       'bootstrap_method':'Class-stratified test-set resampling; same draw shared by all windows; percentile CI conditional on fitted models.',
       'bootstrap_seed':20260922,'same_patients_labels_order_across_windows':True,'metric_formula_verified_against_sklearn':True}
    for split in ['internal','external']:
        yy=None;ids=None;probs=[]
        for h in HOURS:
            folder=root/f'window_{h:02d}h'
            p=pd.read_pickle(folder/f'private_predictions_{split}.pkl');d=pd.read_pickle(folder/'data'/f'{split}.pkl')
            assert np.array_equal(p.stay_id,d.stay_id) and np.array_equal(p.groupHPD,d.groupHPD)
            if yy is None:yy=p.groupHPD.to_numpy(dtype=int);ids=p.stay_id.to_numpy()
            assert np.array_equal(yy,p.groupHPD) and np.array_equal(ids,p.stay_id)
            pr=p[['p_DR','p_RR','p_PW']].to_numpy();probs.append(pr)
            row=tab[(tab.split==split)&(tab.hours==h)].iloc[0]
            m=metrics(yy,pr);assert all(abs(m[k]-row[k])<1e-10 for k in m)
            assert np.allclose(four_auc(yy,pr),[m[k] for k in NAMES],atol=1e-12)
            lock=json.loads((folder/'model_locked.json').read_text());assert (folder/f'private_predictions_{split}.pkl').stat().st_mtime>=lock['time']
            verify['recomputed_evaluations'].append({'split':split,'hours':h,'n':len(yy)})
        groups=[np.flatnonzero(yy==k) for k in [1,2,3]];rng=np.random.default_rng(20260922)
        boot=np.empty((draws,6,4))
        for b in range(draws):
            ix=np.concatenate([rng.choice(g,size=len(g),replace=True) for g in groups])
            y=yy[ix]
            for j in range(6):boot[b,j]=four_auc(y,probs[j][ix])
        points=np.asarray([four_auc(yy,p) for p in probs]);lo,hi=np.quantile(boot,[.025,.975],axis=0)
        for j,h in enumerate(HOURS):
            for k,metric in enumerate(NAMES):cis.append({'split':split,'hours':h,'metric':metric,'n':len(yy),'auc':points[j,k],'lower95':lo[j,k],'upper95':hi[j,k]})
        for j,h in enumerate(HOURS):
            if h==24:continue
            delta=boot[:,j]-boot[:,3];dl,du=np.quantile(delta,[.025,.975],axis=0)
            for k,metric in enumerate(NAMES):changes.append({'split':split,'hours':h,'reference_hours':24,'metric':metric,'delta_auc':points[j,k]-points[3,k],'lower95':dl[k],'upper95':du[k]})
        print('Verified and bootstrapped',split,flush=True)
    pd.DataFrame(cis).to_csv(root/'auc_confidence_intervals.csv',index=False)
    pd.DataFrame(changes).to_csv(root/'paired_changes_vs24h.csv',index=False)
    (root/'evaluation_verification.json').write_text(json.dumps(verify,indent=2))
    (root/'ANALYSIS_COMPLETE.json').write_text(json.dumps({'time':time.time(),'n_windows':6,'bootstrap_draws':draws},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--draws',type=int,default=1000)
    a=p.parse_args();analyze(a.root,a.draws)
