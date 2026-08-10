# Deep K=2 versus K=3 stability experiment

**Status: DEEP_K3_INITIALIZATION_SENSITIVITY_REMAINS**

All 18 prespecified deep fits were retained: three cohorts × K=2/3 × three
independent initializations. Each fit used burn 50, keep 2,000, and thin 50.
Numeric labels were aligned before exact-agreement calculations. Conventional
untreated multi-chain R-hat was not calculated because mixture-label switching
precludes direct pooling.

## Cohort-level K=3 reproducibility

| cohort   |   k3_fits |   k3_fits_passing |   median_agreement_vs_archived |   minimum_ari_vs_archived |   minimum_cluster_prevalence |   maximum_high_lag1_fraction |   minimum_pairwise_agreement |   minimum_pairwise_ari | k3_cross_seed_status               |   K2_selections |   K3_selections |
|:---------|----------:|------------------:|-------------------------------:|--------------------------:|-----------------------------:|-----------------------------:|-----------------------------:|-----------------------:|:-----------------------------------|----------------:|----------------:|
| aumc     |         3 |                 3 |                         0.9977 |                    0.7764 |                       0.0852 |                       0      |                       0.8887 |                 0.774  | PASS_ALL_3_INITIALIZATIONS         |               0 |               3 |
| eicu     |         3 |                 1 |                         0.7706 |                    0.4952 |                       0.0028 |                       0.6667 |                       0.7685 |                 0.4903 | CAUTION_INITIALIZATION_SENSITIVITY |               2 |               1 |
| mimic    |         3 |                 3 |                         0.9924 |                    0.6115 |                       0.0755 |                       0      |                       0.8597 |                 0.6126 | PASS_ALL_3_INITIALIZATIONS         |               0 |               3 |

## K=3 seed versus archived assignments

| cohort   |     seed |   patients | label_mapping            |   exact_agreement_vs_archived |   ari_vs_archived |   nmi_vs_archived |   minimum_cluster_prevalence |
|:---------|---------:|-----------:|:-------------------------|------------------------------:|------------------:|------------------:|-----------------------------:|
| mimic    | 20260805 |       4713 | {"1": 1, "2": 2, "3": 3} |                        0.8593 |            0.6115 |            0.5143 |                       0.0995 |
| mimic    | 20260806 |       4713 | {"1": 1, "2": 2, "3": 3} |                        0.9987 |            0.9972 |            0.9891 |                       0.076  |
| mimic    | 20260807 |       4713 | {"1": 1, "2": 2, "3": 3} |                        0.9924 |            0.9787 |            0.9493 |                       0.0755 |
| eicu     | 20260805 |       1417 | {"1": 1, "2": 2, "3": 3} |                        0.9965 |            0.9901 |            0.9744 |                       0.0882 |
| eicu     | 20260806 |       1417 | {"1": 3, "2": 2, "3": 1} |                        0.7706 |            0.4952 |            0.432  |                       0.0028 |
| eicu     | 20260807 |       1417 | {"1": 3, "2": 2, "3": 1} |                        0.7706 |            0.4952 |            0.432  |                       0.0028 |
| aumc     | 20260805 |       2183 | {"1": 1, "2": 2, "3": 3} |                        0.9977 |            0.9928 |            0.9828 |                       0.0999 |
| aumc     | 20260806 |       2183 | {"1": 1, "2": 2, "3": 3} |                        0.8891 |            0.7764 |            0.6193 |                       0.0852 |
| aumc     | 20260807 |       2183 | {"1": 1, "2": 2, "3": 3} |                        0.9986 |            0.9951 |            0.9894 |                       0.1003 |

## Pairwise K=3 seed stability

| cohort   |   left_seed |   right_seed |   patients | label_mapping            |   exact_agreement |    ari |    nmi |
|:---------|------------:|-------------:|-----------:|:-------------------------|------------------:|-------:|-------:|
| mimic    |    20260805 |     20260806 |       4713 | {"1": 1, "2": 2, "3": 3} |            0.8597 | 0.6126 | 0.5155 |
| mimic    |    20260805 |     20260807 |       4713 | {"1": 1, "2": 2, "3": 3} |            0.8648 | 0.6216 | 0.5303 |
| mimic    |    20260806 |     20260807 |       4713 | {"1": 1, "2": 2, "3": 3} |            0.9928 | 0.9803 | 0.9515 |
| eicu     |    20260805 |     20260806 |       1417 | {"1": 3, "2": 2, "3": 1} |            0.7685 | 0.4903 | 0.4294 |
| eicu     |    20260805 |     20260807 |       1417 | {"1": 3, "2": 2, "3": 1} |            0.7685 | 0.4903 | 0.4294 |
| eicu     |    20260806 |     20260807 |       1417 | {"1": 1, "2": 2, "3": 3} |            1      | 1      | 1      |
| aumc     |    20260805 |     20260806 |       2183 | {"1": 1, "2": 2, "3": 3} |            0.8887 | 0.774  | 0.6155 |
| aumc     |    20260805 |     20260807 |       2183 | {"1": 1, "2": 2, "3": 3} |            0.9991 | 0.9977 | 0.9922 |
| aumc     |    20260806 |     20260807 |       2183 | {"1": 1, "2": 2, "3": 3} |            0.8887 | 0.7751 | 0.6168 |

Every prespecified initialization, including poorly mixing or degenerate fits, is
reported. The experiment evaluates whether a three-pattern representation can be
reproduced; it does not establish K=3 as the unique true taxonomy. The deep-fit
composite score uses the same normalization scope as the fresh grid: deviance and
lag-1 failure are min-max scaled within each cohort and initialization before K=2
and K=3 are compared.
