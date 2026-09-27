"""Fixed-parameter, original-label 6-36h classifiers using pre-MI trajectories."""
import argparse, hashlib, json, time
from pathlib import Path
import numpy as np
import pandas as pd
from autogluon.tabular import TabularPredictor
from prepare import BASE,SPECS,EXPECTED_LABEL_HASH,sha
from train import metrics

BASELINE = {
 'mimic':('00.data_mimic/disease_definition/AKI/df_base_crea.csv','stay_id','baseline_Scr',1.0),
 'eicu':('00.data_eicu/disease_definition/AKI/df_base_crea.csv','stay_id','baseline_creatinine',1.0),
 'aumcdb':('00.data_aumc/disease_definition/AKI/baseline_creatinine.csv','admissionid','baseline_creatinine',.01131),
}
HP={'n_estimators':1000,'max_depth':3,'learning_rate':.03,'min_child_weight':1,
    'reg_lambda':1,'subsample':.85,'colsample_bytree':.85}

def features(tr,ids,endpoint):
    tr=tr[tr.time.isin(range(1,endpoint+1)) & tr.stay_id.isin(ids)].copy()
    assert not tr.duplicated(['stay_id','time']).any()
    g=tr.groupby('stay_id',sort=True);cols={}
    for f in BASE:
        for stat in ['max','min','mean']:cols[f'{f}_{stat}']=getattr(g[f],stat)()
    for i in range(1,endpoint+1):
        for j in range(i+1,endpoint+1):
            d=tr[tr.time.eq(j)].set_index('stay_id')[BASE]-tr[tr.time.eq(i)].set_index('stay_id')[BASE]
            for f in BASE:cols[f'{f}_diff_{i}_{j}']=d[f]
    return pd.DataFrame(cols).reindex(ids).replace([np.inf,-np.inf],np.nan)

def prepare(source,root):
    manifest={'sources':{},'cohorts':{},'windows':{},'temporal_boundary_checks':[],
      'preprocessing':'Read archived forward-filled data BEFORE cohort-wide multiple imputation. Forward carry may include pre-onset observations already available; no backward fill or pooled MICE. Remaining NaNs enter model pipeline.',
      'eligibility':'Fixed original study population and archived partitions at every horizon; no new requirement to survive or remain in ICU to 36h.',
      'parameters':HP,'seed':20260921,'hours':[6,12,18,24,30,36],
      'design':'Fixed hyperparameters from previous 24h training CV; five-fold bagging, no ensemble/stacking, no window chosen on external outcomes.'}
    cohorts={}
    for co,(folder,oldtraj,labfile,n) in SPECS.items():
        lp=source/folder/labfile;fp=source/folder/'df_im_By_ffill.csv'
        assert sha(lp)==EXPECTED_LABEL_HASH[co]
        lab=pd.read_csv(lp,usecols=['stay_id','groupHPD']).astype('int64')
        assert lab.groupby('stay_id').groupHPD.nunique().max()==1
        lab=lab.drop_duplicates('stay_id').set_index('stay_id');assert len(lab)==n
        tr=pd.read_csv(fp);tr['stay_id']=tr.stay_id.astype('int64')
        if co=='mimic':
            tr=tr.rename(columns={c:c[:-5] for c in tr if c.endswith('_mean')})
        if 'bilirubin' not in tr and 'bilirubin_total' in tr:tr['bilirubin']=tr.bilirubin_total
        bp,idcol,bcol,scale=BASELINE[co];bp=source/bp
        baseline=pd.read_csv(bp,usecols=[idcol,bcol]).rename(columns={idcol:'stay_id',bcol:'baseline'})
        baseline['stay_id']=baseline.stay_id.astype('int64');baseline['baseline']*=scale
        tr=tr.merge(baseline,on='stay_id',how='left',validate='many_to_one')
        assert tr.loc[tr.stay_id.isin(lab.index),'baseline'].notna().all()
        tr['crea_divide_basecrea']=(tr.creatinine/tr.baseline).round(2)
        for f in BASE:tr[f]=pd.to_numeric(tr[f],errors='coerce') if f in tr else np.nan
        tr=tr[['stay_id','time']+BASE]
        # Validate the baseline source against ALL original archived ratios.
        old=pd.read_csv(source/folder/oldtraj).merge(baseline,on='stay_id',validate='many_to_one')
        assert old.stay_id.nunique()==n
        assert np.allclose((old.creatinine/old.baseline).round(2),old.crea_divide_basecrea,atol=1e-8)
        cohorts[co]=(tr,lab)
        for p in [lp,fp,bp,source/folder/'step1_数据准备.ipynb']:manifest['sources'][str(p)]=sha(p)
        manifest['cohorts'][co]={'n':n,'baseline_ratio_matches_archive':True,
                                 'class_counts':lab.groupHPD.value_counts().sort_index().to_dict()}
    oldinputs=source/'07.autogluon/01.model/Result-a1234_selfv2_MimiceICU_AUMC_CorrMICfilt/input'
    splits={}
    for split,file in [('train','train_set.csv'),('internal','test_set1.csv'),('external','test_set2.csv')]:
        p=oldinputs/file;manifest['sources'][str(p)]=sha(p)
        splits[split]=pd.read_csv(p,usecols=['stay_id','groupHPD']).astype('int64').set_index('stay_id')
    assert not any(set(splits[a].index)&set(splits[b].index) for a,b in [('train','internal'),('train','external'),('internal','external')])
    for h in manifest['hours']:
        endpoint=h//6;frames=[];obs={}
        for co,(tr,lab) in cohorts.items():
            x=features(tr,lab.index,endpoint)
            # Future perturbation must not affect any feature at this horizon.
            perturbed=tr.copy();perturbed.loc[perturbed.time>endpoint,BASE]=987654.321
            pd.testing.assert_frame_equal(x,features(perturbed,lab.index,endpoint))
            manifest['temporal_boundary_checks'].append({'cohort':co,'hours':h,'future_perturbation_invariant':True})
            x['groupHPD']=lab.groupHPD;x['dataset']=co;frames.append(x)
            counts=tr[tr.time.between(1,endpoint)].groupby('stay_id').size().reindex(lab.index,fill_value=0)
            obs[co]={'patients_with_no_window_rows':int((counts==0).sum()),
                     'patients_with_endpoint_row':int(tr[tr.time.eq(endpoint)].stay_id.isin(lab.index).sum()),
                     'predictor_missing_fraction':float(x.drop(columns=['groupHPD','dataset']).isna().to_numpy().mean())}
        df=pd.concat(frames);assert df.index.is_unique
        folder=root/f'window_{h:02d}h';(folder/'data').mkdir(parents=True,exist_ok=False)
        featurecols=[c for c in df if c not in ['groupHPD','dataset']]
        assert len(featurecols)==27*(3+endpoint*(endpoint-1)//2)
        ms={'n_raw_features':len(featurecols),'splits':{},'observation_coverage':obs}
        for s,ref in splits.items():
            part=df.loc[ref.index].copy()
            assert (part.groupHPD==ref.groupHPD).all()
            assert (part.dataset.eq('aumcdb')).all() if s=='external' else (~part.dataset.eq('aumcdb')).all()
            part.reset_index().to_pickle(folder/'data'/f'{s}.pkl')
            ms['splits'][s]={'n':len(part),'class_counts':part.groupHPD.value_counts().sort_index().to_dict(),
                            'all_predictors_missing_n':int(part[featurecols].isna().all(axis=1).sum()),
                            'missing_cell_fraction':float(part[featurecols].isna().to_numpy().mean())}
        manifest['windows'][str(h)]=ms
    (root/'input_manifest.json').write_text(json.dumps(manifest,indent=2))
    print('PREPARED all 6 windows; original labels, split and temporal boundary checks passed',flush=True)

def fit(root):
    rows=[];verification={};fold_reference=None
    for h in [6,12,18,24,30,36]:
        folder=root/f'window_{h:02d}h'
        assert not (folder/'predictor').exists()
        tr=pd.read_pickle(folder/'data/train.pkl')
        predictor=TabularPredictor(label='groupHPD',problem_type='multiclass',eval_metric='roc_auc_ovo_macro',
          path=str(folder/'predictor'),verbosity=2,learner_kwargs={'random_state':20260921})
        print(f'START {h}h',flush=True)
        predictor.fit(tr.drop(columns=['stay_id','dataset']),hyperparameters={'XGB':HP},
          num_bag_folds=5,num_bag_sets=1,num_stack_levels=0,num_cpus=16,num_gpus=0,
          ag_args_ensemble={'fold_fitting_strategy':'sequential_local'},fit_weighted_ensemble=False,
          calibrate=False,calibrate_decision_threshold=False,time_limit=1200)
        lb=predictor.leaderboard(silent=True);assert len(lb)==1
        lb.to_csv(folder/'training_cv.csv',index=False);model=lb.iloc[0].model
        (folder/'model_locked.json').write_text(json.dumps({'model':model,'time':time.time(),'parameters':HP},indent=2))
        bag=predictor._trainer.load_model(model);assert len(bag.models)==5
        # Actual train/validation fold membership must be identical across windows.
        cv=list(bag._cv_splitters[0].split(tr.drop(columns=['stay_id','dataset','groupHPD']),tr.groupHPD))
        digest=hashlib.sha256(b''.join(np.asarray(i,dtype='int64').tobytes()+np.asarray(j,dtype='int64').tobytes() for i,j in cv)).hexdigest()
        if fold_reference is None:fold_reference=digest
        assert digest==fold_reference
        oof=predictor.predict_proba_oof(model=model).loc[:,[1,2,3]].to_numpy()
        assert abs(metrics(tr.groupHPD.to_numpy(),oof)['auc_ovo_macro']-float(lb.iloc[0].score_val))<1e-10
        details=[]
        for childname in bag.models:
            child=bag.load_child(childname);details.append({'fold':childname,'trained_parameters':child.params_trained})
        (folder/'fit_details.json').write_text(json.dumps({'folds':details,'fold_hash':digest,'model_parameters':predictor.model_hyperparameters(model)},indent=2))
        for s in ['internal','external']:
            test=pd.read_pickle(folder/'data'/f'{s}.pkl');y=test.groupHPD.to_numpy(dtype=int)
            proba=predictor.predict_proba(test.drop(columns=['stay_id','dataset','groupHPD']),model=model).loc[:,[1,2,3]].to_numpy()
            assert np.isfinite(proba).all() and np.allclose(proba.sum(axis=1),1,atol=1e-5)
            pred=pd.DataFrame({'stay_id':test.stay_id,'groupHPD':y,'p_DR':proba[:,0],'p_RR':proba[:,1],'p_PW':proba[:,2]})
            pred.to_pickle(folder/f'private_predictions_{s}.pkl')
            row={'hours':h,'split':s,'n':len(test),'n_raw_features':len(test.columns)-3,'model':model,**metrics(y,proba)}
            rows.append(row);print(json.dumps(row),flush=True)
        verification[str(h)]={'fold_count':len(bag.models),'fold_hash':digest,'oof_auc_recomputed':True}
        pd.DataFrame(rows).to_csv(root/'metrics.csv',index=False)
        (root/'fit_verification.json').write_text(json.dumps(verification,indent=2))
    (root/'TRAINING_COMPLETE.json').write_text(json.dumps({'time':time.time(),'windows':6,'evaluations':len(rows)},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--root',type=Path,required=True)
    a=p.parse_args();a.root.mkdir(exist_ok=False);prepare(a.source,a.root);fit(a.root)
