# Deep multi-seed urine-output documentation sensitivity

**Status: UO_MULTI_SEED_SENSITIVITY_CAUTIONS_REMAIN**

Both reviewer-requested eICU documentation scenarios were run with three deep
independent initializations. Every seed is retained. A scenario passes only when
all three fits satisfy the prespecified mixing, agreement, ARI, and minimum-cluster
prevalence thresholds.

| scenario           |   seeds_completed |   passing_seeds |   patients |   retained_fraction |   median_exact_agreement |   minimum_exact_agreement |   median_ari |   minimum_ari |   median_nmi |   maximum_uncertain_fraction |   minimum_cluster_prevalence |   median_dr_fraction |   median_rr_fraction |   median_pw_fraction |   maximum_high_lag1_fraction | scenario_status                                     |
|:-------------------|------------------:|----------------:|-----------:|--------------------:|-------------------------:|--------------------------:|-------------:|--------------:|-------------:|-----------------------------:|-----------------------------:|---------------------:|---------------------:|---------------------:|-----------------------------:|:----------------------------------------------------|
| documented_windows |                 3 |               3 |       1043 |              0.7361 |                   0.8351 |                    0.8313 |       0.5229 |        0.5128 |       0.4383 |                       0.2972 |                       0.0671 |               0.1774 |               0.7555 |               0.0671 |                       0      | PASS_ALL_3_INITIALIZATIONS                          |
| high_coverage      |                 3 |               0 |        415 |              0.2929 |                   0.8145 |                    0.8048 |       0.5414 |        0.5156 |       0.4014 |                       0.0337 |                       0      |               0.2819 |               0.7181 |               0      |                       0.4167 | CAUTION_DOCUMENTATION_SELECTION_OR_MODE_SENSITIVITY |

## Seed-level evidence

| scenario           |     seed |   patients |   retained_fraction | label_mapping            |   exact_agreement |    ari |    nmi |   uncertain_fraction |   minimum_cluster_prevalence |   dr_fraction |   rr_fraction |   pw_fraction |   high_absolute_lag1_fraction | fit_status   |
|:-------------------|---------:|-----------:|--------------------:|:-------------------------|------------------:|-------:|-------:|---------------------:|-----------------------------:|--------------:|--------------:|--------------:|------------------------------:|:-------------|
| documented_windows | 20260805 |       1043 |              0.7361 | {"1": 1, "2": 2, "3": 3} |            0.836  | 0.5251 | 0.4437 |               0.2637 |                       0.0671 |        0.1774 |        0.7555 |        0.0671 |                        0      | PASS         |
| documented_windows | 20260806 |       1043 |              0.7361 | {"1": 1, "2": 2, "3": 3} |            0.8313 | 0.5128 | 0.4344 |               0.1467 |                       0.0671 |        0.1745 |        0.7584 |        0.0671 |                        0      | PASS         |
| documented_windows | 20260807 |       1043 |              0.7361 | {"1": 1, "2": 2, "3": 3} |            0.8351 | 0.5229 | 0.4383 |               0.2972 |                       0.0681 |        0.1774 |        0.7546 |        0.0681 |                        0      | PASS         |
| high_coverage      | 20260805 |        415 |              0.2929 | {"1": 3, "2": 2, "3": 1} |            0.8048 | 0.5156 | 0.3839 |               0.0337 |                       0      |        0.2602 |        0.7398 |        0      |                        0.4167 | CAUTION      |
| high_coverage      | 20260806 |        415 |              0.2929 | {"1": 3, "2": 2, "3": 1} |            0.8145 | 0.5414 | 0.4014 |               0.0313 |                       0      |        0.2916 |        0.7084 |        0      |                        0.4167 | CAUTION      |
| high_coverage      | 20260807 |        415 |              0.2929 | {"1": 3, "2": 2, "3": 1} |            0.8193 | 0.5572 | 0.4205 |               0.0265 |                       0      |        0.2819 |        0.7181 |        0      |                        0.4167 | CAUTION      |

Restriction to documented or high-coverage windows changes the target population;
agreement in these selected subsets does not establish robustness for all 1,417
original eICU patients.
