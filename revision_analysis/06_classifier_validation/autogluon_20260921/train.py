"""Prespecified CPU parameter comparison; holdouts are opened only after locking selection."""
import argparse
import hashlib
import importlib.metadata as md
import json
import os
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
from autogluon.tabular import TabularPredictor
from sklearn.metrics import roc_auc_score,accuracy_score,balanced_accuracy_score,f1_score,log_loss,confusion_matrix

def config():
    # Explicit grid, same five-fold splits for all candidates; no external-data tuning.
    hp={'XGB':[], 'GBM':[], 'CAT':[], 'XT':[]}
    for depth,rate,child,reg in [(3,.03,1,1),(3,.08,1,1),(4,.04,3,3),(5,.03,5,5),(6,.03,8,10),(2,.05,5,5)]:
        hp['XGB'].append(dict(n_estimators=1000,max_depth=depth,learning_rate=rate,min_child_weight=child,
                              reg_lambda=reg,subsample=.85,colsample_bytree=.85,
                              ag_args={'name_suffix':f'_d{depth}_lr{rate}_mc{child}','priority':100}))
    for leaves,rate,child,reg in [(15,.03,30,1),(31,.03,30,3),(15,.06,60,5),(63,.025,40,5)]:
        hp['GBM'].append(dict(num_boost_round=1500,num_leaves=leaves,learning_rate=rate,
                             min_data_in_leaf=child,lambda_l2=reg,feature_fraction=.85,
                             ag_args={'name_suffix':f'_l{leaves}_lr{rate}','priority':90}))
    for depth,reg in [(4,3),(6,5),(7,10)]:
        hp['CAT'].append(dict(iterations=1200,depth=depth,learning_rate=.04,l2_leaf_reg=reg,
                             ag_args={'name_suffix':f'_d{depth}','priority':80}))
    for leaf in [1,5]:
        hp['XT'].append(dict(n_estimators=500,min_samples_leaf=leaf,max_features=.8,
                            ag_args_ensemble={'use_child_oof':False},
                            ag_args={'name_suffix':f'_leaf{leaf}','priority':70}))
    return {'hyperparameters':hp,'eval_metric':'roc_auc_ovo_macro','seed':20260921,
            'num_bag_folds':5,'num_bag_sets':1,'num_stack_levels':0,'num_cpus':16,
            'time_limit_seconds':1800,'selection':'5-fold out-of-fold macro OVO AUC on original training patients',
            'locked_test_models':['best_XGBoost_by_training_OOF','best_overall_by_training_OOF'],
            'scope':'Original manuscript labels/split; 0-24h; all 243 raw predictors; no new clustering.'}

def metrics(y,p):
    pred=np.argmax(p,axis=1)+1
    d={'auc_ovo_macro':float(roc_auc_score(y,p,multi_class='ovo',average='macro',labels=[1,2,3])),
       'auc_ovr_macro':float(roc_auc_score(y,p,multi_class='ovr',average='macro',labels=[1,2,3])),
       'accuracy':float(accuracy_score(y,pred)),'balanced_accuracy':float(balanced_accuracy_score(y,pred)),
       'f1_macro':float(f1_score(y,pred,average='macro')),'log_loss':float(log_loss(y,p,labels=[1,2,3])),
       'brier_multiclass':float(np.square(p-np.eye(3)[y-1]).sum(axis=1).mean())}
    for a,b,n in [(1,2,'DR_RR'),(1,3,'DR_PW'),(2,3,'RR_PW')]:
        mask=np.isin(y,[a,b]);pa=p[mask,a-1];pb=p[mask,b-1]
        d['auc_'+n]=float(roc_auc_score(y[mask]==a,pa/np.maximum(pa+pb,1e-15)))
    for k,n in [(1,'DR'),(2,'RR'),(3,'PW')]:
        d['auc_ovr_'+n]=float(roc_auc_score(y==k,p[:,k-1]))
        d['sensitivity_'+n]=float((pred[y==k]==k).mean())
    return d

def run(root):
    cfg=config();out=root/'results';out.mkdir(exist_ok=False)
    (out/'plan.json').write_text(json.dumps(cfg,indent=2))
    versions={p:md.version(p) for p in ['autogluon.tabular','autogluon.core','pandas','numpy','scikit-learn','xgboost-cpu','lightgbm','catboost']}
    (out/'environment.json').write_text(json.dumps({'versions':versions,'python':platform.python_version()},indent=2))
    train=pd.read_pickle(root/'data/train.pkl').drop(columns=['stay_id','dataset'])
    np.random.seed(cfg['seed'])
    predictor=TabularPredictor(label='groupHPD',problem_type='multiclass',eval_metric=cfg['eval_metric'],path=str(out/'predictor'),verbosity=2,
                              learner_kwargs={'random_state':cfg['seed']})
    predictor.fit(train_data=train,hyperparameters=cfg['hyperparameters'],num_bag_folds=5,num_bag_sets=1,
                  num_stack_levels=0,num_cpus=16,num_gpus=0,fit_strategy='sequential',
                  ag_args_ensemble={'fold_fitting_strategy':'sequential_local'},
                  time_limit=cfg['time_limit_seconds'],fit_weighted_ensemble=True,
                  calibrate=False,calibrate_decision_threshold=False,raise_on_no_models_fitted=True)
    lb=predictor.leaderboard(silent=True).sort_values('score_val',ascending=False)
    lb.to_csv(out/'training_cv_leaderboard.csv',index=False)
    xgb=lb[lb.model.str.startswith('XGBoost')]
    assert len(xgb),'No fitted XGBoost model'
    locked={'best_XGBoost':xgb.iloc[0].model,'best_overall':lb.iloc[0].model}
    (out/'locked_selection.json').write_text(json.dumps({'selected':locked,'timestamp':time.time(),'basis':cfg['selection']},indent=2))
    (out/'model_hyperparameters.json').write_text(json.dumps({m:predictor.model_hyperparameters(m) for m in set(locked.values())},indent=2,default=str))
    # Internal/external labels have not been passed to fit/leaderboard/model selection.
    rows=[];cms={};cal=[]
    for split in ['internal','external']:
        test=pd.read_pickle(root/f'data/{split}.pkl');y=test.groupHPD.to_numpy(dtype=int)
        features=test.drop(columns=['stay_id','dataset','groupHPD'])
        for role,model in locked.items():
            probs=predictor.predict_proba(features,model=model).loc[:,[1,2,3]].to_numpy()
            assert np.isfinite(probs).all() and np.allclose(probs.sum(axis=1),1,atol=1e-5)
            rows.append({'role':role,'model':model,'split':split,'n':len(test),**metrics(y,probs)})
            # Private predictions support reproducible evaluation and CIs, never public delivery.
            pd.DataFrame({'stay_id':test.stay_id,'groupHPD':y,'p_DR':probs[:,0],'p_RR':probs[:,1],'p_PW':probs[:,2]}).to_pickle(out/f'private_predictions_{split}_{role}.pkl')
            cms[split+'_'+role]=confusion_matrix(y,np.argmax(probs,axis=1)+1,labels=[1,2,3]).tolist()
            for k in range(3):
                b=np.minimum((probs[:,k]*10).astype(int),9)
                for j in range(10):
                    mask=b==j
                    if mask.any(): cal.append({'split':split,'role':role,'class':[1,2,3][k],'bin':j,'n':int(mask.sum()),'mean_prediction':float(probs[mask,k].mean()),'observed_fraction':float((y[mask]==k+1).mean())})
    pd.DataFrame(rows).to_csv(out/'test_metrics.csv',index=False)
    pd.DataFrame(cal).to_csv(out/'calibration_bins.csv',index=False)
    (out/'confusion_matrices.json').write_text(json.dumps(cms,indent=2))
    predictor.model_failures().to_csv(out/'model_failures.csv',index=False)
    (out/'COMPLETE.json').write_text(json.dumps({'completed_at':time.time(),'selected':locked,'n_evaluations':len(rows)},indent=2))
    print(pd.DataFrame(rows).to_string(index=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);a=p.parse_args();run(a.root)
