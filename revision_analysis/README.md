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
| `09_revision_documents/` | Fail-closed document, numerical-consistency, and reviewer-traceability QA | All comments |

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
6. `launch_deep_refits.sh` runs the prespecified 42-fit confirmation grid in a separate environment: 18 K=2/K=3 primary fits, 18 reruns of all screening-depth CAUTION scenarios, and 6 eICU urine-documentation fits.
7. `summarize_deep_k_stability.py`, `summarize_deep_robustness.py`, and `summarize_uo_multiseed_sensitivity.py` require every prespecified initialization and keep failed or degenerate solutions in the denominator.
8. `plot_k_stability_figure.py` regenerates Figure S2 from traceable K=2–5 screening evidence and the deep K=2/K=3 experiment; `verify_revision_documents.py` also requires Table S1 to match its source CSV cell for cell and verifies all seven embedded revised figures (S2, S10, S11a/b, and S12–S14) against their released source assets by SHA-256 or pixel identity.
9. `audit_numerical_consistency.py` maps aggregate sources to manuscript, supplement, response-letter, table, and figure locations; `audit_reviewer_traceability.py` checks that every editor/reviewer item retains its original comment, response, result, limitation, change location, and evidence link.
10. `publication_figure_style.py` fixes journal-width typography and semantic colors; `audit_revision_figures.py` fails closed unless Figures S2 and S10–S14 have all four required export formats and a completed visual-review record bound to the current PNG SHA-256, reviewer, and timestamp.

Patient-level refit inputs and assignments are written only under
`02_revision_outputs/intermediate/`, which is ignored by Git. Aggregate diagnostics
and reports are eligible for the public evidence bundle after a privacy-column audit.

`07_tables_figures/publish_aggregate_evidence.py` publishes only an explicit whitelist
of identifier-free reports, tables, and figures. The publisher stops if any required
file is missing or if a CSV/JSON contains a prohibited patient identifier field.
The whitelist includes the Table S3 audit, deep clustering summaries, adjusted
outcomes, classifier validation, exploratory diuretic evidence, and final document-QA
status, numerical-consistency status, reviewer-traceability status, and figure-QA
status plus identifier-free publication assets, but never
patient-level refit assignments, predictions, manuscripts, or detailed text-extraction
logs.
