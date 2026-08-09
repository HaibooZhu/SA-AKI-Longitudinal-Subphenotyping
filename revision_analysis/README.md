# Revision Analysis Workstreams

_Version-controlled code and documentation for all analyses added during the JTIM major revision._

---

## 📋 Workstream map

| Directory | Purpose | Primary review mapping |
| --- | --- | --- |
| `00_data_lineage/` | Authoritative cohort manifest, artifact hashes and cross-file lineage checks | E.1, R1-M3 |
| `00_config/` | Paths, environment, seeds and shared definitions | All comments |
| `01_data_audit/` | Table S3 and end-to-end data-integrity audit | E.1, R1-M3 |
| `02_missingness_sensitivity/` | Urine output, baseline SCr, RRT and preprocessing sensitivity | E.2, R1-M4, R1-m1 |
| `03_cluster_robustness/` | Cluster selection, stability and assignment uncertainty | E.2, R1-m4 |
| `04_independent_outcomes/` | Adjusted outcomes independent of cluster-defining features | E.3 |
| `05_diuretic_exploratory/` | Reframed post-exposure association analyses | E.4, R1-M1 |
| `06_classifier_validation/` | Balanced performance, calibration and leakage checks | E.5, R1-M2 |
| `07_tables_figures/` | Scripted regeneration of revision tables and figures | All analytical comments |
| `08_literature_update/` | Reproducible PubMed search and literature-screening records | R1-m2 |

## ✍️ File conventions

- Prefix scripts with the mapped comment ID when practical, for example `E1_R1M3_audit_table_s3.py`
- Store configuration rather than absolute paths inside analysis scripts
- Send generated outputs to `../02_revision_outputs`
- Send execution logs to `../03_logs`
- Record input snapshot, configuration, seed, output path and verification status

## Editor Concern #2 implementation

The current revision-only robustness path is intentionally separate from frozen code:

1. `build_eicu_uo_scenarios.py` reconstructs both sum and mean from raw eICU urine-output records and stops unless one rule passes prespecified match, MAE, median-error, and dominance thresholds.
2. `build_cross_cohort_scenarios.py` generates z-score, median/IQR, documented-RRT exclusion, complete 30-window follow-up, and limited-forward-fill/complete-row inputs for all three cohorts.
3. `run_mixak_refit.R` uses a dynamic cohort-specific response list, native 0–1 posterior probabilities, and one random-intercept flag per response.
4. `summarize_fresh_k_grid.py` requires the complete 3-cohort × K=2–5 × 3-seed grid before reporting a result.
5. `summarize_robustness_refits.py` aligns numeric labels before calculating agreement, ARI, NMI, cluster prevalence, and posterior uncertainty.

Patient-level refit inputs and assignments are written only under
`02_revision_outputs/intermediate/`, which is ignored by Git. Aggregate diagnostics
and reports are eligible for the public evidence bundle after a privacy-column audit.

`07_tables_figures/publish_aggregate_evidence.py` publishes only an explicit whitelist
of identifier-free reports, tables, and figures. The publisher stops if any required
file is missing or if a CSV/JSON contains a prohibited patient identifier field.
