# Input data contract

The public repository does not distribute any clinical data. This document records the minimum shapes expected by the revision scripts so authorized researchers can map their local, derived tables without changing source code.

## General rules

- `stay_id` is the internal episode key used to join derived tables. It must remain inside the authorized computing environment.
- One `stay_id` must refer to one ICU episode within a cohort.
- Time-window variables use the study's 6-hour index relative to SA-AKI onset unless a script states otherwise.
- Phenotype labels use `groupHPD`: `1 = DR`, `2 = RR`, `3 = PW`, and `4 = posterior-uncertain` where applicable.
- Units must be harmonized before combining cohorts. The table-regeneration script records the explicit FiO2, bilirubin, and hematocrit rules it applies.

## Revision inputs

| Workstream | Required input | Minimum contract |
| --- | --- | --- |
| Table audit | longitudinal source matrix | `dataset`, `stay_id`, `time`, `groupHPD`, and the features declared in `audit_table_s3.py` |
| Table audit | generated cohort tables | one two-row-header CSV per cohort using the filenames declared in `COHORTS` |
| Workbook extraction | supplementary `.xlsx` files | each audited sheet stores the relevant table in `B2:Z39` |
| Urine sensitivity | event table | `stay_id`, `charttime` in minutes, `urineoutput` |
| Urine sensitivity | onset table | `stay_id`, `saki_onset` in hours |
| Urine sensitivity | clustering input | `stay_id`, `time`, four clustering features, and archived labels |
| Cluster comparison | archived eICU longitudinal table | `stay_id`, `groupHPD` plus four clustering features for plotting |
| Adjusted outcomes | cohort risk-factor tables | files named `df_<cohort>_c<1-3>_riskfactor.csv`; model variables are declared in `FORMULA` |
| Diuretic audit | archived cohort folders | `01.mimic/` and `02.aumc/`, each with the three filenames used by `load_files()` |
| Classifier validation | archived model directory | `input/train_set.csv`, `input/test_set1.csv`, `input/test_set2.csv`, plus archived result summaries used by `original_tables()` |
| Table regeneration | longitudinal source matrix | `dataset`, `stay_id`, `time`, `groupHPD`, and all variables declared in `VARIABLES` |

## Local directory recommendation

```text
authorized_sa_aki_data/
├── source_matrices/
├── cohort_exports/
├── archived_model_inputs/
├── supplementary_workbooks/
└── derived_revision_inputs/
```

Keep this directory outside the Git repository. Pass its files to the scripts through command-line arguments. Store generated outputs under the ignored `results/` directory or another restricted location.

## Pre-run checks

1. Confirm that all source tables come from the frozen analysis snapshot intended for the revision.
2. Verify unique keys and join cardinality before merging.
3. Confirm units and allowed ranges cohort by cohort.
4. Record row counts, patient counts, checksums, code commit, parameters, and random seed.
5. Inspect every candidate Git commit for patient-level files or local credentials.
