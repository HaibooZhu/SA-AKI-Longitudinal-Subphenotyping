"""Restore submitted Table 1 domains from original upstream files, on authoritative labels.

Only aggregate summaries leave this script. The ambiguous legacy AKD summary is
represented by explicitly defined recorded-stage rows in Table S6, not inferred.
"""
from pathlib import Path
import hashlib,json
import pandas as pd
import sys
# (a temporary local dependency directory was added to sys.path here; install the requirements instead)
import numpy as np
from scipy import stats
pvalues=[]
R=Path(__file__).resolve().parents[2];F=R/'00_frozen_inputs/data_snapshot/remote_project_snapshot';O=R/'02_revision_outputs/reports/W25_table1_comparisons_20260908';O.mkdir(exist_ok=True)
meta=[];checks=[]
def read(p):
 p=F/p;meta.append({'source':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()});return pd.read_csv(p)
def unique(d,cols,ids):
 d=d[d.stay_id.isin(ids)][['stay_id']+cols].drop_duplicates()
 assert not d.groupby('stay_id')[cols].nunique(dropna=True).gt(1).any().any(),cols
 return d.groupby('stay_id',as_index=False).first()
def val(s,kind,x=1):
 if kind=='recorded':
  n=len(s);count=int(s.eq(1).sum());return f'{count}/{n} ({100*count/n:.1f})'
 s=s.dropna();n=len(s)
 if not n:return 'Unavailable'
 if kind=='n':return f'{int((s==x).sum())}/{n} ({(s==x).mean()*100:.1f})'
 return f'{s.median():.1f} [{s.quantile(.25):.1f}, {s.quantile(.75):.1f}] (n={n})'
configs=[('mimic','01.MIMICIV_SAKI_trajCluster','sk_first_and_max_stage.csv','lifesupport.csv','sk_event_time.csv'),('aumc','02.AUMCdb_SAKI_trajCluster','aumc_first_and_max_stage.csv','aumcdb_lifesupport.csv','aumcdb_sk_event_time.csv'),('eicu','03.eICU_SAKI_trajCluster','eicu_sk_first_and_max_stage.csv','eicu_lifesupport.csv','eicu_saki_event_time.csv')]
for co,folder,stagefile,supportfile,eventfile in configs:
 d=read(folder+'/sk_survival.csv')[['stay_id','groupHPD','mortality_28d','mortality_7d']];ids=set(d.stay_id);N=len(d);assert len(ids)==N
 def add(s,cols):
  global d
  d=d.merge(unique(s,cols,ids),on='stay_id',how='left',validate='one_to_one');assert len(d)==N
 dm=read(f'00.data_{co}/feature_data/df_{co}_basicinfo.csv');cols=['age','gender','icu_stay_days','icu_expire_flag','outtime']
 if co!='aumc':cols+=['hospital_expire_flag']
 if co=='mimic':cols+=['discharge_location']
 add(dm,cols);d['male']=d.gender.map({'M':1,'F':0})
 bp='baseline_creatinine.csv' if co=='aumc' else 'df_base_crea.csv'
 bs=read(f'00.data_{co}/disease_definition/AKI/'+bp).rename(columns={'admissionid':'stay_id','baseline_creatinine':'baseline_Scr'})
 if co=='aumc':bs.baseline_Scr*=.01131
 add(bs,['baseline_Scr'])
 st=read(f'00.data_{co}/disease_definition/AKI/'+stagefile);add(st,['first_aki_onset','first_aki_stage','aki_endtime_H7D'])
 if co=='mimic':
  d['icu_after']=(pd.to_datetime(d.outtime)-pd.to_datetime(d.first_aki_onset)).dt.total_seconds()/86400
  d['aki_course']=(pd.to_datetime(d.aki_endtime_H7D)-pd.to_datetime(d.first_aki_onset)).dt.total_seconds()/86400
 else:
  d['icu_after']=(d.outtime-d.first_aki_onset)/24;d['aki_course']=(d.aki_endtime_H7D-d.first_aki_onset)/24
 add(read(f'00.data_{co}/disease_definition/AKI/df_peak_discharge.csv'),['Peak_Scr','Discharge_Scr'])
 add(read(f'00.data_{co}/disease_definition/AKI/AKI_diagnose_criteria.csv'),['AKI_criteria'])
 ev=read(f'00.data_{co}/disease_definition/AKI/'+eventfile)
 if co=='mimic':ev['onset_interval']=(pd.to_datetime(ev.aki_onset)-pd.to_datetime(ev.sepsis_onset)).dt.total_seconds()/3600
 else:ev['onset_interval']=ev.aki_onset-ev.sepsis_onset
 add(ev,['onset_interval']);add(read(f'00.data_{co}/treatment/'+supportfile),['is_rrt','is_mv','is_vaso'])
 if co=='aumc':d['hospital_expire_flag']=float('nan')
 groups={'Overall':d,'RR':d[d.groupHPD==2],'DR':d[d.groupHPD==1],'PW':d[d.groupHPD==3]}
 rows=[{'Characteristic':'Patients, N',**{g:str(len(x)) for g,x in groups.items()}}]
 specs=[('Male sex, n/N (%)','male','n',1),('Age, years, median [IQR]','age','med',0),('ICU stay, days, median [IQR]','icu_stay_days','med',0),('ICU stay after SA-AKI onset, days, median [IQR]','icu_after','med',0),('Recorded AKI episode through day 7, days, median [IQR]','aki_course','med',0),('ICU mortality, n/N (%)','icu_expire_flag','n',1),('Hospital mortality, n/N (%)','hospital_expire_flag','n',1),('28-day mortality, n/N (%)','mortality_28d','n',1),('7-day mortality, n/N (%)','mortality_7d','n',1),('Hours from sepsis to AKI onset, median [IQR]','onset_interval','med',0),('Baseline SCr, mg/dL, median [IQR]','baseline_Scr','med',0),('Peak SCr, mg/dL, median [IQR]','Peak_Scr','med',0),('Discharge SCr, mg/dL, median [IQR]','Discharge_Scr','med',0)]
 specs += [(f'AKI criteria: {label}, n/N (%)','AKI_criteria','n',v) for label,v in [('both SCr and urine output','Both'),('SCr only','Scr'),('urine output only','Uo')]]
 specs += [(f'Onset AKI stage {i}, n/N (%)','first_aki_stage','n',i) for i in [1,2,3]]
 specs += [(f'{label}, n/N (%)',col,'recorded',1) for label,col in [('Recorded RRT','is_rrt'),('Recorded mechanical ventilation','is_mv'),('Recorded vasopressor use','is_vaso')]]
 for label,col,kind,v in specs:rows.append({'Characteristic':label,**{g:val(x[col],kind,v) for g,x in groups.items()}})
 assert d.AKI_criteria.notna().all() and set(d.AKI_criteria)=={'Both','Scr','Uo'}
 assert d.first_aki_stage.notna().all() and d.first_aki_stage.isin([1,2,3]).all()
 # Independently reconcile fields shared with the earlier audited W23 summaries.
 core=pd.read_csv(R/f'02_revision_outputs/reports/W23_revision_table1_core_20260908/Table_1_core_{co}.csv').set_index('Characteristic')
 for row in rows:
  if row['Characteristic'] in core.index:
   for g in groups:assert row[g]==str(core.loc[row['Characteristic'],g]),(co,row['Characteristic'],g)
 # Original Table 1 uses non-parametric continuous comparisons and categorical tests.
 def compare(col,kind):
  dd=d[['groupHPD',col]].copy()
  if kind=='recorded':dd[col]=dd[col].eq(1).astype(int)
  dd=dd.dropna(); ng=dd.groupby('groupHPD').size().to_dict()
  if len(ng)<3:return dict(pvalue=None,test='Unavailable',n_by_group=ng)
  if kind=='med':
   arrays=[dd.loc[dd.groupHPD.eq(g),col].to_numpy() for g in [2,1,3]]
   result=stats.kruskal(*arrays);return dict(pvalue=float(result.pvalue),test='Kruskal-Wallis',n_by_group=ng)
  ct=pd.crosstab(dd[col],dd.groupHPD).reindex(columns=[2,1,3],fill_value=0)
  if len(ct)<2:return dict(pvalue=None,test='No variation',n_by_group=ng)
  chi=stats.chi2_contingency(ct,correction=False)
  if chi.expected_freq.min()<5:
   res=stats.fisher_exact(ct.to_numpy(),method=stats.MonteCarloMethod(n_resamples=99999,rng=np.random.default_rng(20260908)))
   return dict(pvalue=float(res.pvalue),test='Fisher-Freeman-Halton Monte Carlo 99999',n_by_group=ng,expected_min=float(chi.expected_freq.min()))
  return dict(pvalue=float(chi.pvalue),test='Pearson chi-square',n_by_group=ng,expected_min=float(chi.expected_freq.min()))
 for label,col,kind,v in specs:
  if any(x['cohort']==co and x['variable']==col for x in pvalues):continue
  result=compare(col,kind);pvalues.append(dict(cohort=co,variable=col,**result))
 pd.DataFrame(rows).to_csv(O/f'Table_1_full_{co}.csv',index=False)
 if co=='mimic':
  remap={'rehabilitation':'Rehabilitation','Another hospital':'Other hospital','death':'Death + Hospice','Hospice':'Death + Hospice'}
  d['discharge_category']=d.discharge_location.replace(remap)
  groups={'Overall':d,'RR':d[d.groupHPD==2],'DR':d[d.groupHPD==1],'PW':d[d.groupHPD==3]}
  cats=['Home','Home health care services','Rehabilitation','Nursing Facility','Other hospital','Death + Hospice','Other']
  assert set(d.discharge_category.dropna()).issubset(cats)
  rows2=[{'Hospital discharge location':c,**{g:val(x.discharge_category,'n',c) for g,x in groups.items()}} for c in cats]
  pd.DataFrame(rows2).to_csv(O/'Table_1_mimic_discharge.csv',index=False)
  assert sum((d.discharge_category==x).sum() for x in cats)==d.discharge_category.notna().sum()
 if co=='mimic':pvalues.append(dict(cohort=co,variable='discharge_category',**compare('discharge_category','n')))
 if co=='eicu':
  patients=read('00.data_eicu/raw/patient.csv')[['patientunitstayid','hospitaldischargelocation']].rename(columns={'patientunitstayid':'stay_id'})
  add(patients,['hospitaldischargelocation'])
  pvalues.append(dict(cohort=co,variable='discharge_category',**compare('hospitaldischargelocation','n')))
 checks.append({'cohort':co,'patients':N,'no_cohort_loss':len(d)==N,'aki_category_exhaustive':True,'onset_stage_exhaustive':True,'shared_W23_fields_match':True,'missing_by_field':{label:int(d[col].isna().sum()) for label,col,kind,v in specs}})
 print(co,N,'all shared W23 cells agree; full rows',len(rows))
(O/'source_manifest.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2));(O/'verification.json').write_text(json.dumps(checks,indent=2))
(O/'README.md').write_text('The original Table 1 domains were reconstructed from the frozen upstream files used by the original cohort-characterization code, joined to the authoritative phenotype cohort. Available-observation denominators are explicit. AKI diagnostic criteria retain the original source labels Both/Scr/Uo. RRT, ventilation and vasopressor rows count positive source indicators among all analyzed cohort members; they describe recorded intervention use, not standardized treatment windows or confirmed non-use when a record is absent. Hospital outcomes remain unavailable in AUMC. The original ambiguous AKD row is replaced by explicitly defined recorded-stage summaries in Table S6; no AKD rate is inferred. The previously corrected eICU discharge categories remain from W22 raw-patient reconstruction. No patient-level files are written.\n')

(O/'table1_comparison_results.json').write_text(json.dumps(pvalues,indent=2))
print('Verified comparisons',len(pvalues),'aggregate only; source input cohort sizes preserved')
