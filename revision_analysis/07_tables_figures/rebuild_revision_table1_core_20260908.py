"""Rebuild a focused Table 1 and AKI-stage table from W5-audited upstream sources.

Avoids mixing legacy characteristics exports with incompatible AKD/AKI-criteria
semantics. Patient data remain private; output consists only of aggregate cells.
"""
from pathlib import Path
import hashlib,json
import pandas as pd
ANA=Path(__file__).resolve().parents[2]
ROOT=ANA/'00_frozen_inputs/data_snapshot/remote_project_snapshot'
OUT=ANA/'02_revision_outputs/reports/W23_revision_table1_core_20260908'
OUT.mkdir(exist_ok=True)
COHORTS={'mimic':('01.MIMICIV_SAKI_trajCluster',4713),'aumc':('02.AUMCdb_SAKI_trajCluster',2183),'eicu':('03.eICU_SAKI_trajCluster',1417)}
sources=[]
def read(p):
    sources.append({'path':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
    return pd.read_csv(p)
def unique(df,cols):
    x=df[['stay_id']+cols].drop_duplicates()
    for c in cols:assert not x.groupby('stay_id')[c].nunique().gt(1).any(),c
    return x.groupby('stay_id',as_index=False).first()
def cell(s,kind,value=1):
    s=s.dropna();n=len(s)
    if not n:return 'Unavailable'
    if kind=='count':return f'{int((s==value).sum())}/{n} ({100*(s==value).mean():.1f})'
    return f'{s.median():.1f} [{s.quantile(.25):.1f}, {s.quantile(.75):.1f}] (n={n})'
for co,(folder,n) in COHORTS.items():
    base=read(ROOT/folder/'sk_survival.csv')[['stay_id','groupHPD','mortality_28d','mortality_7d']]
    assert len(base)==n and not base.stay_id.duplicated().any()
    demos=unique(read(ROOT/f'00.data_{co}/feature_data/df_{co}_basicinfo.csv'),['gender','age'])
    if co=='aumc':
        bp='00.data_aumc/disease_definition/AKI/baseline_creatinine.csv';bc='baseline_creatinine';bid='admissionid';mult=.01131
        sp='00.data_aumc/disease_definition/AKI/aumc_first_and_max_stage.csv'
    elif co=='eicu':
        bp='00.data_eicu/disease_definition/AKI/df_base_crea.csv';bc='baseline_creatinine';bid='stay_id';mult=1
        sp='00.data_eicu/disease_definition/AKI/eicu_sk_first_and_max_stage.csv'
    else:
        bp='00.data_mimic/disease_definition/AKI/df_base_crea.csv';bc='baseline_Scr';bid='stay_id';mult=1
        sp='00.data_mimic/disease_definition/AKI/sk_first_and_max_stage.csv'
    baseline=unique(read(ROOT/bp).rename(columns={bid:'stay_id'}),[bc]).rename(columns={bc:'baseline_Scr'})
    baseline['baseline_Scr']=pd.to_numeric(baseline.baseline_Scr,errors='coerce')*mult
    stcols=['first_aki_stage','max_aki_stage','aki_endstage_H7D']
    stage=unique(read(ROOT/sp),stcols)
    frame=base.merge(demos,on='stay_id',how='left',validate='one_to_one').merge(baseline,on='stay_id',how='left',validate='one_to_one').merge(stage,on='stay_id',how='left',validate='one_to_one')
    frame['male']=frame.gender.map({'M':1,'F':0})
    rows=[]
    groups={'Overall':frame,'RR':frame[frame.groupHPD==2],'DR':frame[frame.groupHPD==1],'PW':frame[frame.groupHPD==3]}
    rows.append({'Characteristic':'Patients, N',**{key:str(len(g)) for key,g in groups.items()}})
    for title,col,kind,value in [('Male sex, n/N (%)','male','count',1),('Age, years, median [IQR]','age','continuous',0),('Baseline SCr, mg/dL, median [IQR]','baseline_Scr','continuous',0)]+[(f'Onset AKI stage {i}, n/N (%)','first_aki_stage','count',i) for i in [1,2,3]]+[('28-day mortality, n/N (%)','mortality_28d','count',1),('7-day mortality, n/N (%)','mortality_7d','count',1)]:
        rows.append({'Characteristic':title,**{key:cell(g[col],kind,value) for key,g in groups.items()}})
    pd.DataFrame(rows).to_csv(OUT/f'Table_1_core_{co}.csv',index=False)
    srows=[]
    for label,col,values in [('At onset','first_aki_stage',[1,2,3]),('At maximum severity','max_aki_stage',[1,2,3]),('At day-7 endpoint','aki_endstage_H7D',[0,1,2,3])]:
        for i in values:srows.append({'Stage assessment':label,'AKI stage':i,**{key:cell(g[col],'count',i) for key,g in groups.items()}})
    pd.DataFrame(srows).to_csv(OUT/f'Table_S6_{co}.csv',index=False)
    print(co,'N',n,'stage nonmissing',frame.first_aki_stage.notna().sum(),'deaths28',int(frame.mortality_28d.sum()))
(OUT/'sources.json').write_text(json.dumps(sources,indent=2))
(OUT/'README.md').write_text('# Focused descriptive tables\n\nTable 1 now focuses on cohort size, age, sex, baseline SCr, onset AKI severity, and independent mortality. Percentages explicitly use nonmissing field denominators. Inputs and phenotype authority follow the W5 audited upstream paths; no legacy characteristics table is promoted solely because its total N matches. Table S6 uses the same source onset/maximum/day-7 stage fields. The day-7 endpoint field is a descriptive recorded renal status, not independent clinical validation or a post-day-7 risk set.\n\nThis replaces inconsistent descriptive table assembly, not trajectory labels or mortality regression results. No new models were fitted. eICU discharge categories remain the separate W22 source-derived panel.\n')
