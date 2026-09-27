"""Audit recovered W21 artifacts without refitting or selecting initializations.

All patient-level assignments and component draws stay in ignored private outputs.
Published summaries are aggregate. Completion, stability and diagnostics are distinct.
"""
from pathlib import Path
import csv, hashlib, itertools, json, os, subprocess
import numpy as np
import pandas as pd

ANA=Path(__file__).resolve().parents[2]
BASE=ANA/'02_revision_outputs/reports/W21_burn_in_remediation'
REC=BASE/'recovered_final_20260908'
OUT=BASE/'verified_summary_20260908';OUT.mkdir(exist_ok=True)
SNAP=ANA/'00_frozen_inputs/data_snapshot/remote_project_snapshot'
MATRICES={'mimic':'01.MIMICIV_SAKI_trajCluster/df_mixAK_fea4_C3.csv','eicu':'03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv','aumc':'02.AUMCdb_SAKI_trajCluster/df_mixAK_fea3_C3_aumc.csv'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
remote='''import hashlib,json,pathlib
r=pathlib.Path.home()/"jtim_w21_burnin_20260823"
files={str(p.relative_to(r)):hashlib.sha256(p.read_bytes()).hexdigest() for p in r.rglob("*") if p.is_file()}
jobs=[x.split("|") for x in (r/"joblist.txt").read_text().splitlines() if x]
inputs={x[2]:hashlib.sha256(pathlib.Path(x[2]).read_bytes()).hexdigest() for x in jobs}
print(json.dumps({"files":files,"inputs":inputs}))
'''
res=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8','-o','StrictHostKeyChecking=yes','-p',os.environ.get('JTIM_SSH_PORT','22'),os.environ['JTIM_SSH_HOST'],'python3 -'],input=remote,text=True,capture_output=True,check=True)
proof=json.loads(res.stdout)
for name,h in proof['files'].items():assert sha(REC/name)==h,name
jobs=[x.split('|') for x in (REC/'joblist.txt').read_text().splitlines() if x]
assert len(jobs)==42
config_hash=sha(REC/'code/mixak_converged_config.json')
expected=[];allrows=[];assign={};chains={}
for co,sc,inp,features,k,seed in jobs:
    k=int(k);seed=int(seed);stem=f'{co}__{sc}__K{k}__seed{seed}';expected.append(stem)
    assert (REC/'logs/job_status'/f'{stem}.exit').read_text().strip()=='0',stem
    p=REC/'outputs'/f'{stem}__diagnostics.csv'
    df=pd.read_csv(p);assert len(df)==1
    row=df.iloc[0].to_dict()
    assert (row['cohort'],row['scenario'],int(row['K']),int(row['seed']))==(co,sc,k,seed)
    assert (int(row['burn']),int(row['keep']),int(row['thin']))==(20000,2000,50)
    assert row['config_sha256']==config_hash
    status=json.loads((REC/'outputs'/f'{stem}__status.json').read_text())
    assert status['config_sha256']==config_hash
    a=pd.read_csv(REC/'outputs'/f'{stem}__assignments.csv')
    assert not a.stay_id.duplicated().any() and len(a)==int(row['patients'])
    assert set(a.group_median).issubset(set(range(1,k+1)))
    ch=pd.read_csv(REC/'outputs'/f'{stem}__mu_chain.csv.gz')
    assert len(ch)==2000 and np.isfinite(ch.to_numpy()).all()
    assign[(co,sc,k,seed)]=a;chains[(co,sc,k,seed)]=ch
    row['minimum_prevalence']=min(a.group_median.value_counts(normalize=True).reindex(range(1,k+1),fill_value=0))
    allrows.append(row)
assert len(set(expected))==42
assert len(list((REC/'outputs').glob('*__diagnostics.csv')))==42

def align(ref,cand,k):
    cross=pd.crosstab(cand,ref).reindex(index=range(1,k+1),columns=range(1,k+1),fill_value=0)
    best=max(itertools.permutations(range(k)),key=lambda p:sum(cross.iloc[i,p[i]] for i in range(k)))
    return {i+1:int(best[i])+1 for i in range(k)}

def adjusted_rand_score(ref,cand):
    tab=pd.crosstab(ref,cand).to_numpy(dtype=float);n=tab.sum()
    choose2=lambda x:x*(x-1)/2
    pairs=choose2(n)
    if pairs==0:return 1.0
    within=choose2(tab).sum();a=choose2(tab.sum(axis=1)).sum();b=choose2(tab.sum(axis=0)).sum()
    expected=a*b/pairs;denom=(a+b)/2-expected
    return 1.0 if denom==0 else (within-expected)/denom

assert adjusted_rand_score(pd.Series([1,1,2,2]),pd.Series([2,2,1,1]))==1.0
assert abs(adjusted_rand_score(pd.Series([1,1,2,2]),pd.Series([1,2,1,2]))+0.5)<1e-12

archives={co:pd.read_csv(SNAP/path,usecols=['stay_id','groupHPD']).drop_duplicates('stay_id') for co,path in MATRICES.items()}
group_summary=[];seed_metrics=[];rhat_rows=[]
diagnostics=pd.DataFrame(allrows)
for (co,sc,k),ds in diagnostics.groupby(['cohort','scenario','K']):
    seeds=sorted(ds.seed.astype(int));assert seeds==[20260805,20260806,20260807]
    reference=assign[(co,sc,k,seeds[0])][['stay_id','group_median']]
    aligned_chains=[];agreements=[];aris=[]
    for seed in seeds:
        a=assign[(co,sc,k,seed)]
        m=a.merge(reference,on='stay_id',suffixes=('_candidate','_reference'),validate='one_to_one')
        assert len(m)==len(a)==len(reference)
        mapping=align(m.group_median_reference,m.group_median_candidate,k)
        ch=chains[(co,sc,k,seed)]
        renames={c:'mu'+str(mapping[int(c.split('.')[0][2:])])+'.'+c.split('.')[1] for c in ch.columns}
        aligned_chains.append(ch.rename(columns=renames)[sorted(ch.columns)].to_numpy())
        if k==3:
            aa=a.merge(archives[co],on='stay_id',validate='one_to_one')
            assert len(aa)==len(a)
            mp=align(aa.groupHPD,aa.group_median,3);al=aa.group_median.map(mp)
            agreement=float((al==aa.groupHPD).mean());ari=float(adjusted_rand_score(aa.groupHPD,al))
            agreements.append(agreement);aris.append(ari)
            seed_metrics.append({'cohort':co,'scenario':sc,'K':k,'seed':seed,'patients':len(aa),'agreement_vs_archived':agreement,'ari_vs_archived':ari,'label_mapping':json.dumps(mp)})
    # Classical split R-hat after assignment-based alignment of component names.
    split=np.stack([half for ch in aligned_chains for half in np.split(ch,2)],axis=0)
    n=split.shape[1];within=split.var(axis=1,ddof=1).mean(axis=0)
    between=n*split.mean(axis=1).var(axis=0,ddof=1)
    rhat=np.sqrt(((n-1)/n*within+between/n)/within)
    assert np.isfinite(rhat).all()
    for col,val in zip(sorted(ch.columns),rhat):rhat_rows.append({'cohort':co,'scenario':sc,'K':k,'component_parameter':col,'classic_split_rhat':float(val)})
    stability=[]
    if k==3:
        for _,r in ds.iterrows():
            sm=next(x for x in seed_metrics if (x['cohort'],x['scenario'],x['seed'])==(co,sc,int(r.seed)))
            stability.append(sm['agreement_vs_archived']>=.75 and sm['ari_vs_archived']>=.5 and r.minimum_prevalence>=.03 and r.high_absolute_lag1_fraction<=0)
    group_summary.append({'cohort':co,'scenario':sc,'K':k,'starts':3,'patients':int(ds.patients.iloc[0]),'passing_legacy_stability_checks':sum(stability) if k==3 else 'NA','median_agreement':float(np.median(agreements)) if agreements else np.nan,'minimum_ari':min(aris) if aris else np.nan,'minimum_prevalence':ds.minimum_prevalence.min(),'maximum_high_lag1_fraction':ds.high_absolute_lag1_fraction.max(),'maximum_mu_geweke_fail_fraction':ds.mu_geweke_fail_fraction.max(),'minimum_mu_ESS':ds.mu_minimum_effective_size.min(),'maximum_classic_split_rhat':float(max(rhat))})
diagnostics.to_csv(OUT/'w21_all_42_diagnostics.csv',index=False)
pd.DataFrame(seed_metrics).to_csv(OUT/'w21_k3_assignment_stability.csv',index=False)
summary=pd.DataFrame(group_summary);summary.to_csv(OUT/'w21_all_scenarios_summary.csv',index=False)
pd.DataFrame(rhat_rows).to_csv(OUT/'w21_classic_split_rhat.csv',index=False)
manifest={'recovered_on':'2026-09-08','server_completion':(REC/'logs/ALL_DONE').read_text(),'verified_files':len(proof['files']),'expected_runs':42,'unique_diagnostics':42,'exit_zero':42,'config_sha256':config_hash,'remote_hashes':proof,'script_sha256':sha(Path(__file__)),'boundary':'Artifact and job completion verified. No formal convergence claim follows from exit status or filenames. Classic split R-hat is not rank-normalized; thresholds and all starts are reported. Source inputs were fingerprinted at recovery, not retrospectively asserted to have runtime fingerprints.'}
(OUT/'recovery_verification.json').write_text(json.dumps(manifest,indent=2))
(OUT/'W21_FINAL_ARTIFACT_AUDIT.md').write_text('# W21 final artifact recovery and aggregate verification\n\nRecovered on 2026-09-08. Server logs record 42/42 successful task exits on 2026-08-24 at 11:07:12 UTC. All copied files match server SHA-256 values at recovery. Each run has assignments, 2,000 retained component-mean draws, diagnostics, and matching configuration hashes. No model was rerun.\n\nThe actual settings were burn=20,000, keep=2,000, thin=50. These refer to thinned blocks, giving 1,100,000 sampler scans. Older explanatory text calling these 20,000 unthinned iterations is not used.\n\nAll 42 initializations are retained. Passing legacy agreement/ARI/prevalence/lag checks does not imply formal convergence. The recovered diagnostics and aligned classical split R-hat are reported separately; no selected subset replaces the original main analysis.\n\n'+'```text\n'+summary.to_string(index=False)+'\n```'+'\n')
print(summary.to_string(index=False))
