"""Run ON THE SERVER, inside the window-prediction AutoGluon 1.6.3 environment, e.g.
    python server_verify_and_shap_24h.py --root ~/jtim_windows_20260921
Step 1 reproduces the saved private_predictions_{internal,external}.pkl of the 24 h window with the saved
predictor (must match to ~1e-6; otherwise the predictor on disk is not the prediction-time model).
Step 2 computes TreeSHAP contributions (xgboost pred_contribs, raw per-class margin) of each of the five fold
models on the internal test set, and checks that the fold boosters reproduce the child predict_proba output.
Aggregate CSV/JSON go to <root>/window_24h/shap/; the per-patient contribution array stays in the same
private folder and must not enter any manuscript deliverable."""
import argparse, json, numpy as np, pandas as pd, xgboost as xgb
from pathlib import Path
from autogluon.tabular import TabularPredictor
p=argparse.ArgumentParser(); p.add_argument('--root',type=Path,required=True); a=p.parse_args()
w=a.root/'window_24h'; out=w/'shap'; out.mkdir(exist_ok=True); model='XGBoost_BAG_L1'
predictor=TabularPredictor.load(str(w/'predictor')); rep={'xgboost':xgb.__version__}
for split in ['internal','external']:
    test=pd.read_pickle(w/'data'/f'{split}.pkl'); saved=pd.read_pickle(w/f'private_predictions_{split}.pkl')
    assert (test.stay_id.to_numpy()==saved.stay_id.to_numpy()).all()
    proba=predictor.predict_proba(test.drop(columns=['stay_id','dataset','groupHPD']),model=model).loc[:,[1,2,3]].to_numpy()
    rep[f'predictor_vs_saved_max_abs_diff_{split}']=float(np.abs(proba-saved[['p_DR','p_RR','p_PW']].to_numpy()).max())
print('step1', rep, flush=True)
bag=predictor._trainer.load_model(model)
test=pd.read_pickle(w/'data/internal.pkl'); X=test.drop(columns=['stay_id','dataset','groupHPD'])
saved=pd.read_pickle(w/'private_predictions_internal.pkl')[['p_DR','p_RR','p_PW']].to_numpy()
contribs=[]; child_probs=[]; raw_probs=[]; feats=None; best=[]
for name in bag.models:
    child=bag.load_child(name); feats=list(child._features_internal or child.features)
    Xp=child.preprocess(X.copy()); Xn=np.asarray(Xp,dtype=np.float32) if not isinstance(Xp,xgb.DMatrix) else None
    bst=child.model.get_booster(); it=(0,int(bst.best_iteration)+1); best.append(it[1]-1)
    dm=xgb.DMatrix(Xn) if Xn is not None else Xp
    raw_probs.append(bst.predict(dm,iteration_range=it)); child_probs.append(child.predict_proba(X.copy()))
    contribs.append(bst.predict(dm,pred_contribs=True,iteration_range=it))
rep['child_predict_proba_avg_vs_saved']=float(np.abs(np.mean(child_probs,0)-saved).max())
rep['raw_booster_avg_vs_saved']=float(np.abs(np.mean(raw_probs,0)-saved).max())
rep['raw_booster_vs_child_predict_proba_max']=float(max(np.abs(r-c).max() for r,c in zip(raw_probs,child_probs)))
rep['best_iterations']=best; rep['n_features']=len(feats)
C=np.stack(contribs)  # folds x n x 3 x (feat+1); last column = bias
np.save(out/'shap_contribs_internal_5fold_PRIVATE.npy',C.astype(np.float32))
mean_abs=np.abs(C[:,:,:,:-1]).mean(axis=1)            # folds x class x feat
imp=pd.DataFrame({'feature':feats,'DR':mean_abs[:,0].mean(0),'RR':mean_abs[:,1].mean(0),'PW':mean_abs[:,2].mean(0)})
imp['total']=imp[['DR','RR','PW']].sum(1); imp['fold_sd_total']=mean_abs.sum(1).std(0)
imp=imp.sort_values('total',ascending=False).reset_index(drop=True); imp['rank']=imp.index+1
imp.to_csv(out/'shap_mean_abs_by_feature_24h.csv',index=False)
rep['shap_scale']='per-class raw margin of each fold model (xgboost pred_contribs), mean |SHAP| averaged over five fold models; not a decomposition of the bagged probability'
json.dump(rep,open(out/'shap_verification.json','w'),indent=1); print(json.dumps(rep,indent=1)); print(imp.head(15).to_string(index=False))
