"""Rebuild 0-24h features on the server, retaining manuscript labels and split.

All patient-level outputs remain in the private server run directory.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

BASE = ['aniongap','bilirubin','po2','mbp','calcium','baseexcess','heart_rate',
        'temperature','hematocrit','wbc','crea_divide_basecrea','bicarbonate',
        'chloride','resp_rate','potassium','pco2','sbp','hemoglobin','lactate',
        'ph','spo2','glucose','urineoutput','sodium','fio2','dbp','creatinine']
SPECS = {
    'mimic': ('01.MIMICIV_SAKI_trajCluster', 'sk_feature_timescale_Fb2.csv', 'df_mixAK_fea4_C3.csv', 4713),
    'eicu': ('03.eICU_SAKI_trajCluster', 'sk_feature_timescale_Fb2_eicu.csv', 'df_mixAK_fea4_C3_eicu.csv', 1417),
    'aumcdb': ('02.AUMCdb_SAKI_trajCluster', 'sk_feature_timescale_Fb2_aumc.csv', 'df_mixAK_fea3_C3_aumc.csv', 2183),
}
EXPECTED_LABEL_HASH = {
    'mimic':'aca43c03eda33ae09354f474a10ecbfd3084351cfb60497b484793fb6ef4314c',
    'eicu':'a79251d11608c48c8c02aaadfe0c16c4a54dbf6ca5efaebd5efa34b550bface4',
    'aumcdb':'7bbc9080a59cda324ed14cd7126871a7c165b34fc0660bc4bd3dcb7508af96af',
}
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''): h.update(chunk)
    return h.hexdigest()

def build(root,out):
    out.mkdir(parents=True,exist_ok=False)
    report={'window_hours':[0,24], 'time_indices':[1,2,3,4], 'phenotypes':{'1':'DR','2':'RR','3':'PW'},
            'features':'27 variables, max/min/mean and six pairwise differences = 243 raw predictors',
            'selection':'No outcome-based preselection before cross-validation; AutoGluon handles training-data preprocessing.',
            'inputs':{},'cohorts':{},'splits':{}}
    frames=[]
    for co,(folder,trajname,labelname,n) in SPECS.items():
        tp,lp=root/folder/trajname,root/folder/labelname
        th,lh=sha(tp),sha(lp)
        assert lh==EXPECTED_LABEL_HASH[co],f'Label source changed: {co}'
        report['inputs'][str(tp)]=th;report['inputs'][str(lp)]=lh
        labels=pd.read_csv(lp,usecols=['stay_id','groupHPD']).astype({'stay_id':'int64','groupHPD':'int64'})
        assert labels.groupby('stay_id').groupHPD.nunique().max()==1
        labels=labels.drop_duplicates('stay_id').set_index('stay_id')
        assert len(labels)==n
        tr=pd.read_csv(tp)
        tr['stay_id']=tr.stay_id.astype('int64')
        tr=tr[tr.time.isin([1,2,3,4]) & tr.stay_id.isin(labels.index)].copy()
        if 'bilirubin' not in tr and 'bilirubin_total' in tr: tr['bilirubin']=tr.bilirubin_total
        assert not tr.duplicated(['stay_id','time']).any()
        for f in BASE: tr[f]=pd.to_numeric(tr[f],errors='coerce') if f in tr else np.nan
        g=tr.groupby('stay_id',sort=True)
        cols={}
        for f in BASE:
            for stat in ['max','min','mean']: cols[f'{f}_{stat}']=getattr(g[f],stat)()
        for i in range(1,5):
            for j in range(i+1,5):
                delta=tr[tr.time.eq(j)].set_index('stay_id')[BASE]-tr[tr.time.eq(i)].set_index('stay_id')[BASE]
                for f in BASE: cols[f'{f}_diff_{i}_{j}']=delta[f]
        x=pd.DataFrame(cols).reindex(labels.index)
        assert len(cols)==243 and set(tr.stay_id)==set(labels.index)
        x['groupHPD']=labels.groupHPD;x['dataset']=co
        x=x.reset_index()
        frames.append(x)
        report['cohorts'][co]={'n':n,'class_counts':labels.groupHPD.value_counts().sort_index().to_dict()}
    allx=pd.concat(frames,ignore_index=True)
    assert not allx.stay_id.duplicated().any(), 'IDs overlap between cohorts'
    allx=allx.set_index('stay_id')
    old=root/'07.autogluon/01.model/Result-a1234_selfv2_MimiceICU_AUMC_CorrMICfilt/input'
    memberships={}
    for name,fn,n in [('train','train_set.csv',4903),('internal','test_set1.csv',1227),('external','test_set2.csv',2183)]:
        p=old/fn;report['inputs'][str(p)]=sha(p)
        prior=pd.read_csv(p,usecols=['stay_id','groupHPD']).astype('int64').set_index('stay_id')
        assert len(prior)==n and prior.index.is_unique
        data=allx.loc[prior.index].copy()
        assert (data.groupHPD==prior.groupHPD).all(), 'Archived classifier differs from original labels'
        assert (data.dataset.eq('aumcdb')).all() if name=='external' else (~data.dataset.eq('aumcdb')).all()
        memberships[name]=set(prior.index)
        data.reset_index().to_pickle(out/f'{name}.pkl')
        report['splits'][name]={'n':n,'class_counts':data.groupHPD.value_counts().sort_index().to_dict(),
                              'sha256':sha(out/f'{name}.pkl')}
    assert not any(memberships[a]&memberships[b] for a,b in [('train','internal'),('train','external'),('internal','external')])
    assert set.union(*memberships.values())==set(allx.index)
    (out/'input_manifest.json').write_text(json.dumps(report,indent=2))
    print(json.dumps({'cohorts':report['cohorts'],'splits':report['splits'],'features':243}),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();build(a.source,a.out)
