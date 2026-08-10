# W1 — Supplementary Table S3 data-integrity audit

## Conclusion

**Audit status: `PASS_WITH_EXPECTED_RENDERING_ERROR`**

Source-to-generated validation passed. The embedded workbook contains exactly the 22 prespecified MIMIC/RR/day-7 rendering mismatches and no others. This audit localizes that specific defect to the source-matrix-to-workbook reporting path; broader analytical provenance is evaluated separately.

## Layer-by-layer checks

| Cohort | Source vs generated CSV mismatches | Generated CSV vs embedded workbook mismatches |
|---|---:|---:|
| mimic | 0 | 22 |
| aumcdb | 0 | 0 |
| eicu | 0 | 0 |

## Detected spreadsheet fill series

| Cohort | Group | Day | Feature range | Length | Value range |
|---|---|---:|---|---:|---|
| mimic | RR | 7 | Fio2–pH | 23 | 46.7739–68.7739 |

## Cross-cohort unit flags

These are audit flags, not automatic corrections. They require confirmation against each database's extraction dictionary.

| Feature | Max/min median ratio | MIMIC median | AUMC median | eICU median |
|---|---:|---:|---:|---:|
| Fio2 | 99.76 | 49.06 | 47.38 | 0.4918 |
| Bilirubin | 17.40 | 2.436 | 19.87 | 1.141 |
| Hematocrit | 14.00 | 29.92 | 2.25 | 31.51 |

## Resolution

The table-generation script constructs a patient-day urine-output sum in `df_fea_add`, but then does not use that object. Table S3 instead summarizes `urineoutput` from the original six-hour rows using an arithmetic mean. The final Methods, caption, and regenerated Table S3 therefore consistently define this quantity as mean six-hour-window urine output, not a patient-day sum. This reporting choice has been resolved and no author decision remains pending.

## Files produced

- `table_s3_embedded_vs_generated.csv`: every embedded cell and its generated-CSV reference.
- `table_s3_source_vs_generated.csv`: every generated cell and its source-matrix reconstruction.
- `table_s3_fill_series.csv`: detected +1 spreadsheet sequences.
- `table_s3_unit_flags.csv`: cross-cohort median-ratio screen.
- `table_s3_missingness.csv`: source-matrix missingness before table aggregation.

Source matrix: `<REVISION_REPOSITORY>/00_frozen_inputs/data_snapshot/remote_project_snapshot/04.other_feature_in_three_dataset/00.data_merge/df_saki_timeseries_feature_all.csv`
