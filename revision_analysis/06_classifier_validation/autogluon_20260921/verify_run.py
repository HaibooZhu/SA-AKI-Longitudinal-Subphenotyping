"""Verify saved-model selection, folds, probabilities and reported metrics."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from autogluon.tabular import TabularPredictor
from sklearn.metrics import roc_auc_score
from train import metrics

def verify(root):
    out=root/'results'
    assert (out/'COMPLETE.json').exists()
    predictor=TabularPredictor.load(str(out/'predictor'))
    lb=pd.read_csv(out/'training_cv_leaderboard.csv')
    selection=json.loads((out/'locked_selection.json').read_text())
    locked=selection['selected'];train=pd.read_pickle(root/'data/train.pkl')
    assert locked['best_overall']==lb.sort_values('score_val',ascending=False).iloc[0].model
    assert locked['best_XGBoost']==lb[lb.model.str.startswith('XGBoost')].sort_values('score_val',ascending=False).iloc[0].model
    report={'selection_before_holdout_evaluation':True,'models':{},'evaluations':{},'metric_sanity_tests':True}
    y=np.array([1,2,3,1,2,3])
    assert metrics(y,np.eye(3)[y-1])['auc_ovo_macro']==1.0
    assert metrics(y,np.full((6,3),1/3))['auc_DR_PW']==0.5
    for model in lb.model:
        if model.startswith('WeightedEnsemble'):continue
        bag=predictor._trainer.load_model(model)
        assert len(bag.models)==5, (model,len(bag.models))
        assert bag._n_repeats_finished==1
        assert not bag.params.get('use_child_oof',False)
        report['models'][model]={'folds':len(bag.models),'completed_repeats':bag._n_repeats_finished}
    for model in set(locked.values()):
        p=predictor.predict_proba_oof(model=model).loc[:,[1,2,3]]
        auc=float(roc_auc_score(train.groupHPD,p,multi_class='ovo',average='macro'))
        saved=float(lb.set_index('model').loc[model,'score_val'])
        assert abs(auc-saved)<1e-10,(model,auc,saved)
        report.setdefault('selected_oof_auc',{})[model]=auc
    table=pd.read_csv(out/'test_metrics.csv')
    for row in table.to_dict(orient='records'):
        key=f"{row['split']}_{row['role']}"
        fp=out/f'private_predictions_{key}.pkl'
        assert fp.stat().st_mtime>=selection['timestamp']
        pred=pd.read_pickle(fp);test=pd.read_pickle(root/f"data/{row['split']}.pkl")
        assert np.array_equal(pred.stay_id,test.stay_id)
        assert np.array_equal(pred.groupHPD,test.groupHPD)
        probs=pred[['p_DR','p_RR','p_PW']].to_numpy()
        m=metrics(pred.groupHPD.to_numpy(dtype=int),probs)
        assert all(abs(m[k]-row[k])<1e-10 for k in m)
        report['evaluations'][key]={'n':len(pred),'metrics_recomputed':True,'labels_match_original':True}
    (out/'verification.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();verify(a.root)
