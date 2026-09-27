"""087 version: same content as build_variable_list_20260923_083.py, except that the transform column is
renamed transform_in_data_preparation and, where the Note records a later conversion for Tables S3-S5 and the
classifier input (FiO2, hematocrit), the transform cell points to the Note. Author request 2026-09-23.
083 note: same content as build_variable_list_20260923_080.py; reader-facing wording replaces internal
audit terms (archived / inspected / existing dictionary / frozen matrix) and three study-list columns are renamed
(trajectory_field, reported_unit, availability). Author request 2026-09-23.
080 note: delivered files no longer name the
reference software (author request 2026-09-23). Output files are Reference_Item_List.csv and
Reference_source_rules.json; columns reference_concept / reference_unit / reference_source_callback.
Internal provenance (repository commit and dictionary hashes) stays in the audit manifest.
Build a source-traced field list and a separate reference crosswalk.

Reads notebook code and CSV headers only. Does not read patient observations,
modify frozen inputs, change manuscript results, or execute extraction code.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import subprocess
from pathlib import Path


CONCEPTS = dict(zip(
    'spo2 fio2 po2 resp_rate sodium chloride potassium calcium urineoutput creatinine crea_divide_basecrea bun heart_rate dbp sbp mbp wbc temperature pco2 baseexcess ph aniongap lactate bicarbonate glucose hematocrit hemoglobin platelets bilirubin inr pt ptt alp ast alt'.split(),
    'spo2 fio2 po2 resp na cl k ca urine crea DERIVED bun hr dbp sbp map wbc temp pco2 be ph anion_gap lact bicar glu hct hgb plt bili MISSING pt ptt alp ast alt'.split(),
))
COHORTS = {
    'miiv': ('MIMIC-IV', 'mimic_iv', '01.MIMICIV_SAKI_trajCluster', '00.data_mimic/feature_data/sk_icu_feature.csv', 'sk_feature_timescale_Fb2.csv'),
    'eicu': ('eICU-CRD', 'eicu_crd', '03.eICU_SAKI_trajCluster', '00.data_eicu/feature_data/eicu_data_merge.csv', 'sk_feature_timescale_Fb2_eicu.csv'),
    'aumc': ('AmsterdamUMCdb', 'amsterdamumcdb', '02.AUMCdb_SAKI_trajCluster', '00.data_aumc/feature_data/aumc_icu_feature.csv', 'sk_feature_timescale_Fb2_aumc.csv'),
}
MERGE = '04.other_feature_in_three_dataset/00.data_merge/Three_dataset_merge_6hwin.ipynb'
FEATURE_MIMIC = '00.data_mimic/feature_data/feature_merge.ipynb'
FEATURE_AUMC = '00.data_aumc/feature_data/generate_aumcdb_icu_feature.ipynb'
REFERENCE_STATUS = 'Reference mapping'
AVAILABILITY = {'available in frozen matrix': 'Available', 'not available': 'Not available'}
SOURCE_NOTICE = ('Reference selector rules adapted from the ricu concept dictionary (GPL-3.0-only; '
                 'https://github.com/eth-mds/ricu); dictionary version identified by the commit below.')


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def header(p):
    with p.open() as f:
        return next(csv.reader(f))


def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def git(p, *args):
    return subprocess.check_output(['git', '-C', str(p), *args], text=True).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--workspace', type=Path, required=True)
    ap.add_argument('--reference-dict-repo', type=Path, required=True)
    ap.add_argument('--reference-dict-ref', required=True)
    ap.add_argument('--concept-dict-path', required=True)   # repository-relative path of concept-dict.json
    ap.add_argument('--data-sources-path', required=True)   # repository-relative path of data-sources.json
    ap.add_argument('--output', type=Path, required=True)
    ap.add_argument('--audit', type=Path, required=True)
    a = ap.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    a.audit.mkdir(parents=True, exist_ok=True)
    exp = a.workspace / '02_revision_workspace/JTIM_revision_experiments_20260805'
    frozen = exp / '00_frozen_inputs/data_snapshot/remote_project_snapshot'
    active = a.workspace / '01_manuscript/14、Journal of Translational Internal Medicine/revision_20260908/13_原稿逐句确认'
    prior = active / 'Supporting_materials/JTIM_Cross_Cohort_Data_Dictionary.csv'
    dictionary = list(csv.DictReader(prior.open(encoding='utf-8-sig')))
    assert set(CONCEPTS) == {r['common_variable'] for r in dictionary}
    ep_bytes = subprocess.check_output(['git','-C',str(a.reference_dict_repo),'show',a.reference_dict_ref+':'+a.concept_dict_path])
    ds_bytes = subprocess.check_output(['git','-C',str(a.reference_dict_repo),'show',a.reference_dict_ref+':'+a.data_sources_path])
    concepts = json.loads(ep_bytes)
    sources = {d['name']: d for d in json.loads(ds_bytes)}
    refsha = git(a.reference_dict_repo, 'rev-parse', a.reference_dict_ref)
    windows = Path(__file__).resolve().parents[1] / '06_classifier_validation/autogluon_20260921'
    tree = ast.parse((windows / 'prepare.py').read_text())
    base = next(ast.literal_eval(n.value) for n in tree.body if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'BASE' for t in n.targets))
    assert len(base) == 27
    hdr = {co: header(frozen / v[3]) for co, v in COHORTS.items()}
    thdr = {co: header(frozen / v[2] / v[4]) for co, v in COHORTS.items()}
    lab_hdr = header(frozen / '00.data_mimic/feature_data/icu_lab_feature_merge.csv')
    vit_hdr = header(frozen / '00.data_mimic/feature_data/icu_vitalsign_feature_merge.csv')
    evidence = {}
    source_hashes = {}

    def cells(rel, numbers):
        p = frozen / rel
        nb = json.loads(p.read_text())
        source_hashes[rel] = sha(p)
        evidence[rel] = {str(i): ''.join(nb['cells'][i]['source']) for i in numbers}

    cells(FEATURE_MIMIC, [2, 3])
    cells(FEATURE_AUMC, [1])
    cells(MERGE, [6, 7, 14, 17])
    for co, v in COHORTS.items():
        cells(v[2] + '/step1_数据准备.ipynb', {'miiv': [18, 22, 24, 26, 28, 29, 31, 34], 'eicu': [16, 20, 22, 25, 27, 28, 30], 'aumc': [18, 20, 22, 25, 27, 28, 30]}[co])

    study, items, rule_export = [], [], {}
    direct_pairs = set()
    for d in dictionary:
        var = d['common_variable']
        concept = CONCEPTS[var]
        for co, (name, prefix, folder, infile, trajfile) in COHORTS.items():
            alias = {'platelets': 'platelet', 'bilirubin': 'bilirubin_total'}.get(var, var) if co == 'miiv' else ('platelet' if co == 'aumc' and var == 'platelets' else var)
            sourcefield = alias if alias in hdr[co] else ''
            trajfield = alias if alias in thdr[co] else (var if var in thdr[co] else '')
            origin = infile
            refs = folder + '/step1_数据准备.ipynb'
            aggregation = 'Arithmetic mean of source records in each 6-hour window'
            if co == 'miiv':
                refs += '#cells22,24,34'
                if alias in vit_hdr:
                    origin = '00.data_mimic/feature_data/icu_vitalsign_feature_merge.csv -> ' + infile
                elif alias in lab_hdr:
                    origin = '00.data_mimic/feature_data/icu_lab_feature_merge.csv -> ' + infile
                refs += '; ' + FEATURE_MIMIC + '#cells2,3'
            else:
                refs += '#cells20,22'
                if co == 'aumc':
                    refs += '; ' + FEATURE_AUMC + '#cell1'
            transform = ''
            notes = ''
            if co == 'miiv' and alias != var:
                transform = f'{alias} renamed to {var}'
                refs += '; ' + MERGE + '#cells14,17'
            if co == 'aumc':
                upstream_alias = {'temperature':'temperture','po2':'PaO2','pco2':'PaCO2','fio2':'FiO2','urineoutput':'urine_output'}.get(var)
                if upstream_alias:
                    transform = f'{upstream_alias} renamed to {var}'
                if var in ['pt','calcium']:
                    transform = 'Prepared value multiplied by 10 in the feature-construction code'
                    # 070 (2026-09-23): notes aligned with Table S4 (063) and the external classifier input (AUMC calcium set to missing).
                    notes = {'calcium': 'Values are consistent with ionized calcium in mmol/L multiplied by 10 and are not comparable with total calcium in MIMIC-IV and eICU-CRD. Table S4 reports the values divided by 10 as ionized calcium; all calcium-derived classifier predictors were set to missing in the external AUMC test set.',
                             'pt': 'Values are consistent with INR multiplied by 10. Table S4 reports the values divided by 10 as INR. Not a classifier input.'}[var]
                if var == 'spo2':
                    transform = 'Values <10 multiplied by 100 before filtering'
                if var == 'bilirubin':
                    origin = 'database/AMUCdb/raw/feature/bilirubin.csv (referenced in the analysis code)'
                    sourcefield = 'value'
                    transform = 'value renamed to bilirubin; measuredat divided by 1000*60*60; 6-hour mean'
                    refs = MERGE + '#cell7'
                    notes = 'Added during the descriptive-table merge; this is separate from the classifier pre-MI input.'
                if var == 'bun':
                    notes = 'Not used for AUMC clustering. Presence of a reference item does not change the study input.'
            if var == 'crea_divide_basecrea':
                origin = folder + '/' + trajfile
                sourcefield = 'creatinine / baseline creatinine'
                trajfield = var
                transform = 'Creatinine divided by the cohort-specific baseline creatinine, rounded to 2 decimals'
                aggregation = 'Ratio derived from the 6-hour creatinine mean'
                refs = folder + '/step1_数据准备.ipynb#cell' + ('31' if co == 'miiv' else '30')
                if co == 'aumc':
                    transform += '; baseline creatinine multiplied by 0.01131 before division'
                notes = 'Derived variable; no independent itemid.'
            if co == 'eicu':
                missing = 'Within-patient forward fill, then multiple imputation (MICE); urineoutput first filled with zero' if var == 'urineoutput' else 'Within-patient forward fill, then multiple imputation (MICE)'
            else:
                missing = 'Within-patient forward fill, then multiple imputation (MICE)'
            if not sourcefield and var != 'crea_divide_basecrea':
                notes = (notes + ' ' if notes else '') + 'Not in the cohort feature file.'
            # Unit conversions later applied to the descriptive Tables S3-S5 (048) and,
            # for FiO2 and hematocrit, to the classifier input before feature derivation (050).
            if var == 'fio2':
                rule = ('source values are fractions and were multiplied by 100; ' if co == 'eicu' else '') + 'values outside 21-100% were set to missing.'
                notes = (notes + ' ' if notes else '') + 'Tables S3-S5 and classifier input: ' + rule
            if var == 'hematocrit':
                notes = (notes + ' ' if notes else '') + ('Tables S3-S5 and classifier input: values of 0.05-0.80 were treated as fractions '
                                                          'and multiplied by 100; values outside 5-80% were set to missing.')
            if var == 'bilirubin' and co == 'aumc':
                notes = notes + ' Tables S3-S5: source values in micromol/L were divided by 17.1 to mg/dL.'
            # 087: the transform column covers the data-preparation code; later conversions are in the Note column.
            if not transform and 'Tables S3-S5 and classifier input:' in notes:
                transform = 'No conversion in data preparation; see Note for Tables S3-S5 and classifier input'
            cluster = var in ['creatinine','urineoutput','crea_divide_basecrea'] or (var == 'bun' and co != 'aumc')
            concept_data = concepts.get(concept, {})
            ss = concept_data.get('sources', {}).get(co, [])
            if ss:
                ref_status = REFERENCE_STATUS
                direct_pairs.add((var, co))
            elif concept == 'DERIVED':
                ref_status = 'Study-derived ratio; no independent itemid'
            elif concept == 'anion_gap':
                ref_status = 'Derived concept in the reference dictionary; not proof of the study aniongap calculation'
            else:
                ref_status = 'No matching selector in the reference dictionary'
            study.append({
                'variable':var, 'display_name':d['display_name'], 'cohort':name,
                'clustering_input':'Yes' if cluster else 'No',
                'classifier_27_variable_schema':'Yes' if var in base else 'No',
                'feature_file_field':sourcefield or 'Not in the cohort feature file',
                'trajectory_field':trajfield or 'Not in the cohort trajectory file',
                'reported_unit':d['common_unit'],
                'availability':AVAILABILITY[d[prefix+'_availability']],
                'input_file_or_upstream_reference':origin,
                'six_hour_aggregation':aggregation,
                'transform_in_data_preparation':transform or 'No conversion in the available analysis code',
                'documented_missingness':missing,
                'study_code_reference':refs,
                'reference_concept':concept if concept not in ['DERIVED','MISSING'] else '',
                'item_mapping_status':ref_status,
                'note':notes,
            })
            if ss:
                for n, s in enumerate(ss):
                    defaults = sources[co]['tables'][s['table']].get('defaults', {})
                    value = s.get('val_var', s.get('value_var', defaults.get('val_var','')))
                    index = s.get('index_var', defaults.get('index_var',''))
                    ruleid = f'{concept}.{co}.{n+1}'
                    rule_export[ruleid] = {'concept':concept,'cohort':co,'source_selector':s,'table_defaults':defaults,
                        'concept_unit':concept_data.get('unit'),'concept_min':concept_data.get('min'),'concept_max':concept_data.get('max'),
                        'concept_callback':concept_data.get('callback'),'provenance':REFERENCE_STATUS}
                    vals = s.get('ids')
                    if vals is not None:
                        selectors = [('itemid' if s.get('sub_var')=='itemid' else 'label', str(v)) for v in (vals if isinstance(vals,list) else [vals])]
                    elif s.get('regex'):
                        selectors = [('regex', s['regex'])]
                    else:
                        selectors = [('column', value)]
                    for kind, val in selectors:
                        items.append({'variable':var,'cohort':name,'reference_concept':concept,'source_table':s['table'],
                            'selector_type':kind,'selector_field':s.get('sub_var','') if kind != 'column' else '',
                            'selector_value':val,'value_field':value,'time_field':index,
                            'reference_unit':'; '.join(concept_data['unit']) if isinstance(concept_data.get('unit'),list) else concept_data.get('unit','') or 'unitless',
                            'reference_source_callback':s.get('callback',''),
                            'source_rule_id':ruleid,'mapping_status':REFERENCE_STATUS,
                            'source_json_pointer':f'concept-dict.json#/{concept}/sources/{co}/{n}'})
            elif concept == 'anion_gap':
                ruleid = f'{concept}.{co}.derived'
                rule_export[ruleid] = {'concept':concept,'cohort':co,'definition':concept_data,'provenance':REFERENCE_STATUS}
                items.append({'variable':var,'cohort':name,'reference_concept':concept,'source_table':'Derived concept',
                    'selector_type':'derived','selector_field':'','selector_value':'See anion_gap callback and dependencies in Reference_source_rules.json',
                    'value_field':'','time_field':'','reference_unit':'mEq/L; mmol/l','reference_source_callback':concept_data.get('callback',''),
                    'source_rule_id':ruleid,'mapping_status':ref_status,'source_json_pointer':'concept-dict.json#/anion_gap'})
            else:
                items.append({'variable':var,'cohort':name,'reference_concept':concept if concept not in ['DERIVED','MISSING'] else '',
                    'source_table':'','selector_type':'derived' if concept=='DERIVED' else 'unmapped',
                    'selector_field':'','selector_value':'creatinine / baseline creatinine' if concept=='DERIVED' else '',
                    'value_field':'','time_field':'','reference_unit':'','reference_source_callback':'','source_rule_id':'',
                    'mapping_status':ref_status,'source_json_pointer':''})

    assert len(study)==105 and len({(r['variable'],r['cohort']) for r in study})==105
    assert sum(r['clustering_input']=='Yes' for r in study)==11
    assert sum(r['classifier_27_variable_schema']=='Yes' for r in study)==81
    assert {r['variable'] for r in items}==set(CONCEPTS)
    assert len({(r['variable'],r['cohort']) for r in items})==105
    # Exact comparison against the source dictionary, including every item ID.
    for var, co in direct_pairs:
        expected=[]
        for s in concepts[CONCEPTS[var]]['sources'][co]:
            ids=s.get('ids')
            if ids is not None:
                expected.extend(str(x) for x in (ids if isinstance(ids,list) else [ids]))
            else: expected.append(s.get('regex',s.get('val_var',s.get('value_var',''))))
        actual=[r['selector_value'] for r in items if r['variable']==var and r['cohort']==COHORTS[co][0]]
        assert actual==expected, (var,co)
    write_csv(a.output/'JTIM_Study_Field_List.csv', study)
    write_csv(a.output/'Reference_Item_List.csv', items)
    (a.output/'Reference_source_rules.json').write_text(json.dumps({'source_notice':SOURCE_NOTICE,'reference_dictionary_commit':easysha,'concept_dict_sha256':hashlib.sha256(ep_bytes).hexdigest(),
        'data_sources_sha256':hashlib.sha256(ds_bytes).hexdigest(),'historical_study_extraction_equivalence':'Not established','rules':rule_export},ensure_ascii=False,indent=2))
    (a.audit/'study_code_excerpts.json').write_text(json.dumps(evidence,ensure_ascii=False,indent=2))
    (a.audit/'workbook_input.json').write_text(json.dumps({'study':study,'items':items},ensure_ascii=False))
    counts={'study_variables':35,'cohort_variable_rows':len(study),'direct_reference_cohort_variable_pairs':len(direct_pairs),
        'item_list_rows':len(items),'numeric_itemid_rows':sum(r['selector_type']=='itemid' for r in items),
        'label_rows':sum(r['selector_type']=='label' for r in items),'column_rows':sum(r['selector_type']=='column' for r in items),
        'regex_rows':sum(r['selector_type']=='regex' for r in items),
        'unmapped':[{k:r[k] for k in ['variable','cohort','mapping_status']} for r in items if r['selector_type']=='unmapped'],
        'no_patient_observations_read':True,'checks':['105 unique study rows','35 variables across all three cohorts','11 clustering inputs','27 classifier variables per cohort','reference selectors match dictionary exactly']}
    (a.audit/'verification.json').write_text(json.dumps(counts,ensure_ascii=False,indent=2))
    manifest={'created_on':'2026-09-23','reference_dict_commit':refsha,'reference_concept_dict_sha256':hashlib.sha256(ep_bytes).hexdigest(),'reference_data_sources_sha256':hashlib.sha256(ds_bytes).hexdigest(),
        'previous_dictionary_sha256':sha(prior),'frozen_notebook_sha256':source_hashes,
        'classifier_code_sha256':{f:sha(windows/f) for f in ['prepare.py','run_windows.py']},
        'counts':counts,'outputs':{p.name:sha(p) for p in a.output.iterdir() if p.is_file()}}
    (a.audit/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    print(json.dumps(counts,ensure_ascii=False))


if __name__=='__main__':
    main()
