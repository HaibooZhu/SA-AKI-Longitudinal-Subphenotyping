#!/usr/bin/env python3
"""Evaluate the locked 24-h XGBoost_BAG_L1 model from the private six-window backup:
pairwise ROC curves (pair-normalized scores), class-specific metrics, one-vs-rest calibration,
and fold-model SHAP (raw margin scale). Reads patient-level files only from the private tmp
extraction; writes aggregate CSV/JSON and figures. No patient identifiers are written."""
import os, json, pickle, sys, hashlib
import numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_fscore_support, accuracy_score
import statsmodels.api as sm
from scipy.special import logit
import xgboost as xgb
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','..'))
PRIV=os.path.join(ROOT,'tmp','jtim_windows_extract','jtim_windows_20260921','window_24h')
OUT=os.path.join(ROOT,'02_revision_outputs','reports','20260921_window_prediction','24h_model_evaluation')
os.makedirs(OUT,exist_ok=True)
CLASSES=[1,2,3]; NAME={1:'DR',2:'RR',3:'PW'}; PAIRS=[(1,2,'DR vs RR'),(1,3,'DR vs PW'),(2,3,'RR vs PW')]
COHORT={'internal':'MIMIC-IV + eICU internal test set','external':'AUMC external test set'}

def load_pred(split):
    df=pickle.load(open(os.path.join(PRIV,f'private_predictions_{split}.pkl'),'rb'))
    y=df['groupHPD'].to_numpy().astype(int); P=df[['p_DR','p_RR','p_PW']].to_numpy(float)
    return df, y, P

def calibration_parameters(yb, prob):
    clipped=np.clip(prob,1e-6,1-1e-6); design=sm.add_constant(logit(clipped))
    fit=sm.GLM(yb, design, family=sm.families.Binomial()).fit(); return float(fit.params[0]), float(fit.params[1])

summary={}; roc_rows=[]; class_rows=[]; pair_rows=[]
for split in ['internal','external']:
    df,y,P=load_pred(split)
    macro=float(roc_auc_score(y,P,multi_class='ovo',average='macro',labels=CLASSES))
    summary[split]={'n':int(len(y)),'macro_ovo_auc':macro}
    for a,b,lab in PAIRS:
        m=np.isin(y,[a,b]); pa=P[m,a-1]; pb=P[m,b-1]; s=pa/np.maximum(pa+pb,1e-15)
        auc=float(roc_auc_score(y[m]==a,s)); fpr,tpr,_=roc_curve(y[m]==a,s)
        pair_rows.append({'split':split,'pair':lab,'n':int(m.sum()),'auc_pair_normalized':auc})
        for f,t in zip(fpr,tpr): roc_rows.append({'split':split,'pair':lab,'fpr':float(f),'tpr':float(t)})
    pred=np.array(CLASSES)[P.argmax(1)]
    pr,rc,f1,sup=precision_recall_fscore_support(y,pred,labels=CLASSES,zero_division=0)
    summary[split]['accuracy']=float(accuracy_score(y,pred))
    for i,c in enumerate(CLASSES):
        yb=(y==c).astype(int); ic,sl=calibration_parameters(yb,P[:,i])
        class_rows.append({'split':split,'class':NAME[c],'n':int(sup[i]),'precision':float(pr[i]),'recall':float(rc[i]),'f1':float(f1[i]),
                           'auc_ovr':float(roc_auc_score(yb,P[:,i])),'calibration_intercept':ic,'calibration_slope':sl})
pairs=pd.DataFrame(pair_rows); classes=pd.DataFrame(class_rows); rocs=pd.DataFrame(roc_rows)
pairs.to_csv(os.path.join(OUT,'pairwise_auc_24h.csv'),index=False); classes.to_csv(os.path.join(OUT,'class_metrics_calibration_24h.csv'),index=False)
rocs.to_csv(os.path.join(OUT,'roc_curves_24h.csv'),index=False)
# cross-check against the recorded metrics.csv of the six-window run
rec=pd.read_csv(os.path.join(ROOT,'02_revision_outputs','reports','20260921_window_prediction','metrics.csv')); rec=rec[rec.hours==24]
check={}
for split in ['internal','external']:
    r=rec[rec.split==split].iloc[0]
    check[split]={'macro_recorded':float(r.auc_ovo_macro),'macro_here':summary[split]['macro_ovo_auc'],
                  'pairs_recorded':[float(r.auc_DR_RR),float(r.auc_DR_PW),float(r.auc_RR_PW)],
                  'pairs_here':pairs[pairs.split==split].auc_pair_normalized.tolist(),
                  'accuracy_recorded':float(r.accuracy),'accuracy_here':summary[split]['accuracy']}
print(json.dumps(check,indent=1))
ok=all(abs(check[s]['macro_recorded']-check[s]['macro_here'])<1e-9 and np.allclose(check[s]['pairs_recorded'],check[s]['pairs_here'],atol=1e-9) for s in check)
print('MATCHES RECORDED METRICS:',ok)
print(pairs.round(3).to_string(index=False)); print(classes.round(3).to_string(index=False))

# ---------------- Figure 5: pairwise ROC curves, internal (A) and external (B)
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,2,figsize=(11,5.2))
colors={'DR vs RR':'#1f77b4','DR vs PW':'#d62728','RR vs PW':'#2ca02c'}
for ax,split,letter in zip(axes,['internal','external'],['A','B']):
    for a,b,lab in PAIRS:
        sub=rocs[(rocs.split==split)&(rocs.pair==lab)]; auc=pairs[(pairs.split==split)&(pairs.pair==lab)].auc_pair_normalized.iloc[0]
        ax.plot(sub.fpr,sub.tpr,color=colors[lab],lw=2,label=f'{lab} (AUC {auc:.3f})')
    ax.plot([0,1],[0,1],'k--',lw=1.2)
    macro=summary[split]['macro_ovo_auc']
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate'); ax.set_xlim(-0.01,1.01); ax.set_ylim(-0.01,1.01)
    ax.set_title(f'{COHORT[split]} (n = {summary[split]["n"]:,})\nmacro one-versus-one AUC {macro:.3f}',fontsize=11)
    ax.legend(loc='lower right',frameon=False,title='Subphenotype pair'); ax.set_aspect('equal')
    ax.text(-0.12,1.06,letter,transform=ax.transAxes,fontsize=18,fontweight='bold')
fig.tight_layout()
for ext in ['pdf','svg','png']: fig.savefig(os.path.join(OUT,f'Figure5_pairwise_ROC_24h.{ext}'),dpi=300 if ext=='png' else None,bbox_inches='tight')
fig.savefig(os.path.join(OUT,'Figure5_pairwise_ROC_24h.tiff'),dpi=600,bbox_inches='tight',pil_kwargs={'compression':'tiff_lzw'})
fig.savefig(os.path.join(OUT,'JTIM_Figure_5.jpg'),dpi=300,bbox_inches='tight')
plt.close(fig)

# ---------------- reproduce saved probabilities from the five fold boosters (column alignment + iteration range)
X_train=pickle.load(open(os.path.join(PRIV,'predictor','utils','data','X.pkl'),'rb'))
cols=list(X_train.columns); print('training feature columns:',len(cols))
feat=pickle.load(open(os.path.join(PRIV,'data','internal.pkl'),'rb'))
extra=[c for c in feat.columns if c not in cols]; print('non-feature columns in internal.pkl:',extra)
dfp,y,P=load_pred('internal')
assert (feat['stay_id'].to_numpy()==dfp['stay_id'].to_numpy()).all(), 'row order differs between feature and prediction files'
Xi=feat[cols].astype(float); dm=xgb.DMatrix(Xi.to_numpy(), feature_names=None)
boosters=[]; best=[]
for k in range(1,6):
    b=xgb.Booster(); b.load_model(os.path.join(PRIV,'predictor','models','XGBoost_BAG_L1',f'S1F{k}','xgb.ubj')); boosters.append(b); best.append(int(b.attr('best_iteration')))
def fold_probs(use_best):
    out=[]
    for b,bi in zip(boosters,best):
        out.append(b.predict(dm, iteration_range=(0,bi+1) if use_best else (0,0)))
    return np.mean(out,axis=0)
res={}
for use_best in [True,False]:
    Q=fold_probs(use_best); res[use_best]=float(np.max(np.abs(Q-P)))
print('max |reproduced - saved| probability: best_iteration+1 ->',res[True],'| all trees ->',res[False])
use_best=res[True]<res[False]; print('use best_iteration range:',use_best,'best_iterations:',best)
assert min(res.values())<1e-5, 'could not reproduce saved probabilities'

# ---------------- SHAP contributions on each fold model's raw margin (per class), internal test set
contribs=[]
for b,bi in zip(boosters,best):
    c=b.predict(dm, pred_contribs=True, iteration_range=(0,bi+1) if use_best else (0,0))  # (n, 3, n_feat+1)
    contribs.append(c[:,:,:-1])  # drop bias term
C=np.stack(contribs)  # folds x n x class x feat
mean_abs_fold=np.abs(C).mean(axis=1)            # folds x class x feat
mean_abs=mean_abs_fold.mean(axis=0)             # class x feat  (fold-averaged mean |SHAP|)
imp=pd.DataFrame({'feature':cols,'DR':mean_abs[0],'RR':mean_abs[1],'PW':mean_abs[2]}); imp['total']=imp[['DR','RR','PW']].sum(1)
imp['fold_sd_total']=np.abs(C).mean(axis=1).sum(axis=1).std(axis=0)
imp=imp.sort_values('total',ascending=False).reset_index(drop=True); imp['rank']=imp.index+1
imp.to_csv(os.path.join(OUT,'shap_mean_abs_by_feature_24h.csv'),index=False)
print(imp.head(15).round(4).to_string(index=False))
# save fold-mean per-patient SHAP for plotting only (no identifiers), plus feature values
Cm=C.mean(axis=0)  # n x class x feat (fold-averaged contributions, margin scale)
np.save(os.path.join(ROOT,'tmp','jtim_windows_extract','shap_foldmean_internal.npy'),Cm)
Xi.to_pickle(os.path.join(ROOT,'tmp','jtim_windows_extract','X_internal_aligned.pkl'))
json.dump({'summary':summary,'check':check,'matches_recorded':bool(ok),'reproduction_max_abs_diff':res,'use_best_iteration':bool(use_best),'best_iterations':best,
           'shap_scale':'per-class raw margin of each fold model (xgboost pred_contribs), averaged over the five fold models',
           'internal_features_sha256':hashlib.sha256(Xi.to_numpy().tobytes()).hexdigest()},open(os.path.join(OUT,'evaluation_24h.json'),'w'),indent=1)
print('done')
