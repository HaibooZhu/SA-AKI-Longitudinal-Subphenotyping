# W5 Independent clinical outcome analysis

## Bottom line

This revision analysis uses authoritative phenotype and survival files for every
cohort and adjusts mortality associations for age, sex, AKI stage at first diagnosis,
baseline serum creatinine, and onset-window nonrenal SOFA. The endpoint is independent
of the kidney variables used to define the seven-day phenotypes. Results support
prognostic association only; they do not establish causal treatment response.

The strict day-7 landmark analysis is a sensitivity analysis. Requiring a complete
trajectory through time 28 preferentially retains longer-observed ICU stays and can
introduce selection/collider bias. Therefore, the all-day-7-survivor analysis is shown
alongside it, and neither replaces the full-cohort mortality analysis.

## Authoritative cohort and lineage checks

Phenotype, 28-day mortality, and 7-day mortality were overwritten from each cohort's
frozen `sk_survival.csv`. Age, sex, baseline serum creatinine, first AKI stage,
and RRT were rebuilt directly from their frozen upstream source files. The legacy
risk-factor exports were used only as comparison targets and did not supply model
values.

| cohort   |   authoritative_n |   legacy_risk_n |   legacy_extra_n |   legacy_group_mismatch_n |   legacy_mortality_mismatch_n |   authoritative_mortality28_missing_n |   authoritative_mortality7_missing_n |
|:---------|------------------:|----------------:|-----------------:|--------------------------:|------------------------------:|--------------------------------------:|-------------------------------------:|
| MIMIC-IV |              4713 |            4713 |                0 |                         0 |                             0 |                                     0 |                                    0 |
| eICU-CRD |              1417 |            1748 |              331 |                       148 |                             0 |                                    17 |                                    7 |
| AUMC     |              2183 |            2183 |                0 |                         0 |                             0 |                                     0 |                                    0 |

Aggregate upstream-versus-legacy covariate verification:

| cohort   | covariate       | upstream_source                                                     |   authoritative_n |   upstream_nonmissing_n |   legacy_nonmissing_n |   comparable_n |   legacy_vs_upstream_mismatch_n | revision_value_source   | missing_value_rule                                              |
|:---------|:----------------|:--------------------------------------------------------------------|------------------:|------------------------:|----------------------:|---------------:|--------------------------------:|:------------------------|:----------------------------------------------------------------|
| MIMIC-IV | gender          | 00.data_mimic/feature_data/df_mimic_basicinfo.csv                   |              4713 |                    4713 |                  4713 |           4713 |                               0 | frozen upstream file    | preserved as missing                                            |
| MIMIC-IV | age             | 00.data_mimic/feature_data/df_mimic_basicinfo.csv                   |              4713 |                    4713 |                  4713 |           4713 |                               0 | frozen upstream file    | preserved as missing                                            |
| MIMIC-IV | baseline_Scr    | 00.data_mimic/disease_definition/AKI/df_base_crea.csv               |              4713 |                    4713 |                  4713 |           4713 |                               0 | frozen upstream file    | preserved as missing                                            |
| MIMIC-IV | first_aki_stage | 00.data_mimic/disease_definition/AKI/sk_first_and_max_stage.csv     |              4713 |                    4713 |                  4713 |           4713 |                               0 | frozen upstream file    | preserved as missing                                            |
| MIMIC-IV | is_rrt          | 00.data_mimic/treatment/lifesupport.csv                             |              4713 |                    4559 |                  4559 |           4559 |                               0 | frozen upstream file    | preserved as missing                                            |
| eICU-CRD | gender          | 00.data_eicu/feature_data/df_eicu_basicinfo.csv                     |              1417 |                    1417 |                  1417 |           1417 |                               0 | frozen upstream file    | preserved as missing                                            |
| eICU-CRD | age             | 00.data_eicu/feature_data/df_eicu_basicinfo.csv                     |              1417 |                    1417 |                  1417 |           1417 |                               0 | frozen upstream file    | preserved as missing                                            |
| eICU-CRD | baseline_Scr    | 00.data_eicu/disease_definition/AKI/df_base_crea.csv                |              1417 |                    1417 |                  1417 |           1417 |                               0 | frozen upstream file    | preserved as missing                                            |
| eICU-CRD | first_aki_stage | 00.data_eicu/disease_definition/AKI/eicu_sk_first_and_max_stage.csv |              1417 |                    1417 |                  1417 |           1417 |                              19 | frozen upstream file    | preserved as missing                                            |
| eICU-CRD | is_rrt          | 00.data_eicu/treatment/eicu_lifesupport.csv                         |              1417 |                    1417 |                  1417 |           1417 |                               0 | frozen upstream file    | left-join absence explicitly set to 0 in archived eICU notebook |
| AUMC     | gender          | 00.data_aumc/feature_data/df_aumc_basicinfo.csv                     |              2183 |                    2134 |                  2134 |           2134 |                               0 | frozen upstream file    | preserved as missing                                            |
| AUMC     | age             | 00.data_aumc/feature_data/df_aumc_basicinfo.csv                     |              2183 |                    2183 |                  2183 |           2183 |                               0 | frozen upstream file    | preserved as missing                                            |
| AUMC     | baseline_Scr    | 00.data_aumc/disease_definition/AKI/baseline_creatinine.csv         |              2183 |                    2183 |                  2183 |           2183 |                               0 | frozen upstream file    | preserved as missing                                            |
| AUMC     | first_aki_stage | 00.data_aumc/disease_definition/AKI/aumc_first_and_max_stage.csv    |              2183 |                    2183 |                  2183 |           2183 |                              34 | frozen upstream file    | preserved as missing                                            |
| AUMC     | is_rrt          | 00.data_aumc/treatment/aumcdb_lifesupport.csv                       |              2183 |                    2117 |                  2117 |           2117 |                               0 | frozen upstream file    | preserved as missing                                            |

The 1,748-row legacy eICU risk-factor export is not an analysis population. It contains
331 patients outside the authoritative 1,417-patient cohort; its legacy phenotype and
mortality columns were ignored.

## Outcome populations

| cohort   | cohort_label   | analysis_population   | population_label                            |   n_selected |   n_with_outcome |   events |   event_risk |   retained_from_authoritative |
|:---------|:---------------|:----------------------|:--------------------------------------------|-------------:|-----------------:|---------:|-------------:|------------------------------:|
| mimic    | MIMIC-IV       | overall               | Full authoritative cohort                   |         4713 |             4713 |      804 |        0.171 |                         1     |
| mimic    | MIMIC-IV       | day7_all_survivors    | Alive at day 7                              |         4277 |             4277 |      368 |        0.086 |                         0.907 |
| mimic    | MIMIC-IV       | day7_full_trajectory  | Alive at day 7 and observed through time 28 |         1237 |             1237 |      269 |        0.217 |                         0.262 |
| eicu     | eICU-CRD       | overall               | Full authoritative cohort                   |         1417 |             1400 |      167 |        0.119 |                         1     |
| eicu     | eICU-CRD       | day7_all_survivors    | Alive at day 7                              |         1309 |             1299 |       66 |        0.051 |                         0.924 |
| eicu     | eICU-CRD       | day7_full_trajectory  | Alive at day 7 and observed through time 28 |          501 |              494 |       50 |        0.101 |                         0.354 |
| aumc     | AUMC           | overall               | Full authoritative cohort                   |         2183 |             2183 |      545 |        0.25  |                         1     |
| aumc     | AUMC           | day7_all_survivors    | Alive at day 7                              |         1892 |             1892 |      254 |        0.134 |                         0.867 |
| aumc     | AUMC           | day7_full_trajectory  | Alive at day 7 and observed through time 28 |         1116 |             1116 |      215 |        0.193 |                         0.511 |

`day7_all_survivors` includes all patients documented alive at day 7 with an observed
28-day outcome. `day7_full_trajectory` additionally requires a kidney trajectory row
at time 28; this is a restrictive sensitivity population, not a less-biased primary
cohort.

## Covariate definition

The primary model is:

`outcome ~ phenotype + age/10 + sex + first AKI stage + log(baseline SCr) + onset nonrenal SOFA`

Onset nonrenal SOFA is the sum of respiratory, coagulation, liver, cardiovascular, and
CNS SOFA components at aligned time 1 (the first six-hour window after SA-AKI onset).
The renal component is deliberately excluded to avoid readjusting for creatinine or
urine output, which define the phenotypes. These archived SOFA components were derived
from rolling 24-hour scores, aligned to six-hour windows, aggregated by the within-window
maximum, and completed within patient by forward/backward filling; residual missing
values were set to zero in AUMC and eICU. This provenance limits causal interpretation.

## Unadjusted 28-day mortality

| cohort_label   | phenotype   |    n |   events |   risk |
|:---------------|:------------|-----:|---------:|-------:|
| MIMIC-IV       | DR          | 1298 |      328 |  0.253 |
| MIMIC-IV       | RR          | 3055 |      282 |  0.092 |
| MIMIC-IV       | PW          |  360 |      194 |  0.539 |
| eICU-CRD       | DR          |  415 |       59 |  0.142 |
| eICU-CRD       | RR          |  862 |       57 |  0.066 |
| eICU-CRD       | PW          |  123 |       51 |  0.415 |
| AUMC           | DR          |  564 |      177 |  0.314 |
| AUMC           | RR          | 1400 |      242 |  0.173 |
| AUMC           | PW          |  219 |      126 |  0.575 |

## Primary adjusted 28-day mortality

Reference phenotype: RR. Cohort-specific binomial GLMs use HC3 robust standard errors.

| cohort_label   | comparison   |   adjusted_or |   ci_low |   ci_high |   p_value |   n_complete |   events |   events_per_parameter |   model_auc | diagnostic_flag   |
|:---------------|:-------------|--------------:|---------:|----------:|----------:|-------------:|---------:|-----------------------:|------------:|:------------------|
| MIMIC-IV       | DR vs RR     |         2.89  |    2.414 |     3.461 |     0     |         4713 |      804 |                 89.333 |       0.738 | PASS              |
| MIMIC-IV       | PW vs RR     |         9.891 |    7.726 |    12.663 |     0     |         4713 |      804 |                 89.333 |       0.738 | PASS              |
| eICU-CRD       | DR vs RR     |         1.932 |    1.293 |     2.889 |     0.001 |         1400 |      167 |                 18.556 |       0.745 | PASS              |
| eICU-CRD       | PW vs RR     |         8.398 |    5.232 |    13.481 |     0     |         1400 |      167 |                 18.556 |       0.745 | PASS              |
| AUMC           | DR vs RR     |         2.032 |    1.595 |     2.589 |     0     |         2134 |      535 |                 59.444 |       0.725 | PASS              |
| AUMC           | PW vs RR     |         5.557 |    4.025 |     7.673 |     0     |         2134 |      535 |                 59.444 |       0.725 | PASS              |

## Minimal-model sensitivity

This model omits nonrenal SOFA and otherwise uses the same baseline covariates.

| cohort_label   | comparison   |   adjusted_or |   ci_low |   ci_high |   p_value |   n_complete |   events |   events_per_parameter |   model_auc | diagnostic_flag   |
|:---------------|:-------------|--------------:|---------:|----------:|----------:|-------------:|---------:|-----------------------:|------------:|:------------------|
| MIMIC-IV       | DR vs RR     |         3.16  |    2.64  |     3.781 |         0 |         4713 |      804 |                100.5   |       0.733 | PASS              |
| MIMIC-IV       | PW vs RR     |        11.681 |    9.146 |    14.919 |         0 |         4713 |      804 |                100.5   |       0.733 | PASS              |
| eICU-CRD       | DR vs RR     |         2.096 |    1.414 |     3.106 |         0 |         1400 |      167 |                 20.875 |       0.722 | PASS              |
| eICU-CRD       | PW vs RR     |        10.779 |    6.821 |    17.035 |         0 |         1400 |      167 |                 20.875 |       0.722 | PASS              |
| AUMC           | DR vs RR     |         2.219 |    1.751 |     2.811 |         0 |         2134 |      535 |                 66.875 |       0.709 | PASS              |
| AUMC           | PW vs RR     |         6.618 |    4.821 |     9.085 |         0 |         2134 |      535 |                 66.875 |       0.709 | PASS              |

## Day-7 landmark analyses: death during days 8–28

Unadjusted event rates:

| cohort_label   | analysis_population   | phenotype   |    n |   events |   risk |
|:---------------|:----------------------|:------------|-----:|---------:|-------:|
| MIMIC-IV       | day7_all_survivors    | DR          | 1143 |      173 |  0.151 |
| MIMIC-IV       | day7_all_survivors    | RR          | 2912 |      139 |  0.048 |
| MIMIC-IV       | day7_all_survivors    | PW          |  222 |       56 |  0.252 |
| MIMIC-IV       | day7_full_trajectory  | DR          |  446 |      120 |  0.269 |
| MIMIC-IV       | day7_full_trajectory  | RR          |  682 |       99 |  0.145 |
| MIMIC-IV       | day7_full_trajectory  | PW          |  109 |       50 |  0.459 |
| eICU-CRD       | day7_all_survivors    | DR          |  388 |       32 |  0.082 |
| eICU-CRD       | day7_all_survivors    | RR          |  828 |       23 |  0.028 |
| eICU-CRD       | day7_all_survivors    | PW          |   83 |       11 |  0.133 |
| eICU-CRD       | day7_full_trajectory  | DR          |  170 |       25 |  0.147 |
| eICU-CRD       | day7_full_trajectory  | RR          |  279 |       18 |  0.065 |
| eICU-CRD       | day7_full_trajectory  | PW          |   45 |        7 |  0.156 |
| AUMC           | day7_all_survivors    | DR          |  494 |      107 |  0.217 |
| AUMC           | day7_all_survivors    | RR          | 1286 |      128 |  0.1   |
| AUMC           | day7_all_survivors    | PW          |  112 |       19 |  0.17  |
| AUMC           | day7_full_trajectory  | DR          |  340 |       94 |  0.276 |
| AUMC           | day7_full_trajectory  | RR          |  698 |      103 |  0.148 |
| AUMC           | day7_full_trajectory  | PW          |   78 |       18 |  0.231 |

Strict full-trajectory sensitivity:

| cohort_label   | comparison   |   adjusted_or |   ci_low |   ci_high |   p_value |   n_complete |   events |   events_per_parameter |   model_auc | diagnostic_flag    |
|:---------------|:-------------|--------------:|---------:|----------:|----------:|-------------:|---------:|-----------------------:|------------:|:-------------------|
| MIMIC-IV       | DR vs RR     |         1.964 |    1.438 |     2.682 |     0     |         1237 |      269 |                 29.889 |       0.692 | PASS               |
| MIMIC-IV       | PW vs RR     |         5.155 |    3.236 |     8.212 |     0     |         1237 |      269 |                 29.889 |       0.692 | PASS               |
| eICU-CRD       | DR vs RR     |         2.205 |    1.129 |     4.305 |     0.021 |          494 |       50 |                  5.556 |       0.686 | CAUTION_LOW_EVENTS |
| eICU-CRD       | PW vs RR     |         2.719 |    1.008 |     7.333 |     0.048 |          494 |       50 |                  5.556 |       0.686 | CAUTION_LOW_EVENTS |
| AUMC           | DR vs RR     |         2.187 |    1.547 |     3.094 |     0     |         1095 |      213 |                 23.667 |       0.678 | PASS               |
| AUMC           | PW vs RR     |         1.8   |    1.012 |     3.201 |     0.046 |         1095 |      213 |                 23.667 |       0.678 | PASS               |

All documented day-7 survivors:

| cohort_label   | comparison   |   adjusted_or |   ci_low |   ci_high |   p_value |   n_complete |   events |   events_per_parameter |   model_auc | diagnostic_flag    |
|:---------------|:-------------|--------------:|---------:|----------:|----------:|-------------:|---------:|-----------------------:|------------:|:-------------------|
| MIMIC-IV       | DR vs RR     |         3.145 |    2.47  |     4.005 |     0     |         4277 |      368 |                 40.889 |       0.724 | PASS               |
| MIMIC-IV       | PW vs RR     |         5.964 |    4.166 |     8.538 |     0     |         4277 |      368 |                 40.889 |       0.724 | PASS               |
| eICU-CRD       | DR vs RR     |         2.623 |    1.454 |     4.731 |     0.001 |         1299 |       66 |                  7.333 |       0.747 | CAUTION_LOW_EVENTS |
| eICU-CRD       | PW vs RR     |         4.459 |    1.988 |     9.997 |     0     |         1299 |       66 |                  7.333 |       0.747 | CAUTION_LOW_EVENTS |
| AUMC           | DR vs RR     |         2.49  |    1.842 |     3.366 |     0     |         1850 |      251 |                 27.889 |       0.697 | PASS               |
| AUMC           | PW vs RR     |         1.703 |    0.988 |     2.934 |     0.055 |         1850 |      251 |                 27.889 |       0.697 | PASS               |

## Covariate-standardized mortality risks

Phenotype was counterfactually set to DR, RR, or PW for every complete-case patient;
predictions were averaged over that cohort's observed covariate distribution. Intervals
are patient-level nonparametric bootstrap percentile intervals.

| cohort_label   | analysis_population   | phenotype   |   standardized_risk |   ci_low |   ci_high |   bootstrap_successes |   bootstrap_failures |
|:---------------|:----------------------|:------------|--------------------:|---------:|----------:|----------------------:|---------------------:|
| MIMIC-IV       | overall               | DR          |               0.235 |    0.209 |     0.257 |                   300 |                    0 |
| MIMIC-IV       | overall               | RR          |               0.098 |    0.088 |     0.108 |                   300 |                    0 |
| MIMIC-IV       | overall               | PW          |               0.503 |    0.455 |     0.55  |                   300 |                    0 |
| MIMIC-IV       | day7_full_trajectory  | DR          |               0.255 |    0.218 |     0.292 |                   300 |                    0 |
| MIMIC-IV       | day7_full_trajectory  | RR          |               0.151 |    0.125 |     0.182 |                   300 |                    0 |
| MIMIC-IV       | day7_full_trajectory  | PW          |               0.463 |    0.37  |     0.56  |                   300 |                    0 |
| eICU-CRD       | overall               | DR          |               0.129 |    0.097 |     0.158 |                   300 |                    0 |
| eICU-CRD       | overall               | RR          |               0.073 |    0.057 |     0.094 |                   300 |                    0 |
| eICU-CRD       | overall               | PW          |               0.373 |    0.291 |     0.465 |                   300 |                    0 |
| eICU-CRD       | day7_full_trajectory  | DR          |               0.136 |    0.084 |     0.18  |                   300 |                    0 |
| eICU-CRD       | day7_full_trajectory  | RR          |               0.068 |    0.038 |     0.101 |                   300 |                    0 |
| eICU-CRD       | day7_full_trajectory  | PW          |               0.162 |    0.058 |     0.298 |                   300 |                    0 |
| AUMC           | overall               | DR          |               0.304 |    0.266 |     0.344 |                   300 |                    0 |
| AUMC           | overall               | RR          |               0.182 |    0.164 |     0.204 |                   300 |                    0 |
| AUMC           | overall               | PW          |               0.525 |    0.465 |     0.602 |                   300 |                    0 |
| AUMC           | day7_full_trajectory  | DR          |               0.273 |    0.227 |     0.326 |                   300 |                    0 |
| AUMC           | day7_full_trajectory  | RR          |               0.15  |    0.123 |     0.176 |                   300 |                    0 |
| AUMC           | day7_full_trajectory  | PW          |               0.238 |    0.152 |     0.331 |                   300 |                    0 |

## Secondary endpoint: in-hospital renal replacement therapy

RRT was not part of the phenotype-defining trajectory, but exact treatment timing is
not available in these derived files. It is therefore secondary and associative.

| cohort_label   | comparison   |   adjusted_or |   ci_low |   ci_high |   p_value |   n_complete |   events |   events_per_parameter |   model_auc | diagnostic_flag    |
|:---------------|:-------------|--------------:|---------:|----------:|----------:|-------------:|---------:|-----------------------:|------------:|:-------------------|
| MIMIC-IV       | DR vs RR     |         7.647 |    4.263 |    13.717 |     0     |         4559 |      168 |                 18.667 |       0.906 | PASS               |
| MIMIC-IV       | PW vs RR     |        46.605 |   26.007 |    83.518 |     0     |         4559 |      168 |                 18.667 |       0.906 | PASS               |
| eICU-CRD       | DR vs RR     |         1.088 |    0.48  |     2.464 |     0.84  |         1417 |       38 |                  4.222 |       0.737 | CAUTION_LOW_EVENTS |
| eICU-CRD       | PW vs RR     |         2.652 |    1.075 |     6.54  |     0.034 |         1417 |       38 |                  4.222 |       0.737 | CAUTION_LOW_EVENTS |
| AUMC           | DR vs RR     |        11.635 |    7.547 |    17.937 |     0     |         2068 |      254 |                 28.222 |       0.863 | PASS               |
| AUMC           | PW vs RR     |        17.466 |   10.801 |    28.244 |     0     |         2068 |      254 |                 28.222 |       0.863 | PASS               |

## Missingness and diagnostics

| cohort   | outcome           | analysis_population   | model_variant               |   n_population |   n_complete |   excluded_missing |   excluded_missing_percent | variables_required                                                                               |
|:---------|:------------------|:----------------------|:----------------------------|---------------:|-------------:|-------------------:|---------------------------:|:-------------------------------------------------------------------------------------------------|
| MIMIC-IV | mortality_28d     | overall               | minimal_baseline            |           4713 |         4713 |                  0 |                      0     | mortality_28d, groupHPD, age10, male, first_aki_stage, log_baseline_scr                          |
| MIMIC-IV | mortality_28d     | overall               | primary_onset_nonrenal_sofa |           4713 |         4713 |                  0 |                      0     | mortality_28d, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa     |
| MIMIC-IV | mortality_day8_28 | day7_all_survivors    | primary_onset_nonrenal_sofa |           4277 |         4277 |                  0 |                      0     | mortality_day8_28, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa |
| MIMIC-IV | mortality_day8_28 | day7_full_trajectory  | primary_onset_nonrenal_sofa |           1237 |         1237 |                  0 |                      0     | mortality_day8_28, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa |
| MIMIC-IV | is_rrt            | overall               | primary_onset_nonrenal_sofa |           4713 |         4559 |                154 |                      3.268 | is_rrt, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa            |
| eICU-CRD | mortality_28d     | overall               | minimal_baseline            |           1417 |         1400 |                 17 |                      1.2   | mortality_28d, groupHPD, age10, male, first_aki_stage, log_baseline_scr                          |
| eICU-CRD | mortality_28d     | overall               | primary_onset_nonrenal_sofa |           1417 |         1400 |                 17 |                      1.2   | mortality_28d, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa     |
| eICU-CRD | mortality_day8_28 | day7_all_survivors    | primary_onset_nonrenal_sofa |           1309 |         1299 |                 10 |                      0.764 | mortality_day8_28, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa |
| eICU-CRD | mortality_day8_28 | day7_full_trajectory  | primary_onset_nonrenal_sofa |            501 |          494 |                  7 |                      1.397 | mortality_day8_28, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa |
| eICU-CRD | is_rrt            | overall               | primary_onset_nonrenal_sofa |           1417 |         1417 |                  0 |                      0     | is_rrt, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa            |
| AUMC     | mortality_28d     | overall               | minimal_baseline            |           2183 |         2134 |                 49 |                      2.245 | mortality_28d, groupHPD, age10, male, first_aki_stage, log_baseline_scr                          |
| AUMC     | mortality_28d     | overall               | primary_onset_nonrenal_sofa |           2183 |         2134 |                 49 |                      2.245 | mortality_28d, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa     |
| AUMC     | mortality_day8_28 | day7_all_survivors    | primary_onset_nonrenal_sofa |           1892 |         1850 |                 42 |                      2.22  | mortality_day8_28, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa |
| AUMC     | mortality_day8_28 | day7_full_trajectory  | primary_onset_nonrenal_sofa |           1116 |         1095 |                 21 |                      1.882 | mortality_day8_28, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa |
| AUMC     | is_rrt            | overall               | primary_onset_nonrenal_sofa |           2183 |         2068 |                115 |                      5.268 | is_rrt, groupHPD, age10, male, first_aki_stage, log_baseline_scr, onset_nonrenal_sofa            |

`CAUTION_LOW_EVENTS` means fewer than 10 outcome events per fitted coefficient. It is a
warning about precision and possible overfitting, not an automatic model failure.

## Reporting decision

1. Use the nonrenal-SOFA-adjusted full-cohort 28-day mortality model as the principal
   independent clinical endpoint.
2. Report both day-7 landmark definitions and explicitly label the full-trajectory
   version as a selection-sensitive analysis.
3. Retain RRT only as a secondary association with a treatment-timing limitation.
4. Do not interpret phenotype coefficients as causal effects or evidence that a
   phenotype-specific treatment improves survival.
