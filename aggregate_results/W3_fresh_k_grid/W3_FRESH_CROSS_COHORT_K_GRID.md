# Fresh cross-cohort K=2–5 screening grid

**Status: COMPLETE_GRID_NONUNANIMOUS_K3**

All 36 prespecified fits were completed: three cohorts × four candidate K values ×
three deterministic seeds. Every fit used the same screening-depth burn, keep, thin,
lag-1 diagnostic, HPD uncertainty rule, and native 0–1 posterior probability scale.

| cohort   |   K |   seeds_completed |    mean_deviance |   sd_deviance |   mean_high_absolute_lag1_fraction |   mean_uncertain_fraction |   median_max_posterior_probability |   median_selection_score |   selections_across_three_seeds |   mean_elapsed_seconds |
|:---------|----:|------------------:|-----------------:|--------------:|-----------------------------------:|--------------------------:|-----------------------------------:|-------------------------:|--------------------------------:|-----------------------:|
| aumc     |   2 |                 3 | 598412           |        1.2003 |                             0      |                    0.0281 |                             0.9973 |                   0.0097 |                               1 |                17.9905 |
| aumc     |   3 |                 3 | 597284           |       73.0614 |                             0.037  |                    0.1345 |                             0.9702 |                   0      |                               2 |                25.4876 |
| aumc     |   4 |                 3 | 712062           |      388.123  |                             0.7083 |                    0.0202 |                             1      |                   1.3903 |                               0 |                34.7017 |
| aumc     |   5 |                 3 | 709706           |     4643.36   |                             0.7111 |                    0.0418 |                             0.9635 |                   1.3793 |                               0 |                37.5665 |
| eicu     |   2 |                 3 | 556129           |       17.5135 |                             0      |                    0.0294 |                             0.9993 |                   0      |                               2 |                14.9933 |
| eicu     |   3 |                 3 | 657124           |    87992.9    |                             0.4306 |                    0.0593 |                             1      |                   1.3017 |                               1 |                26.4801 |
| eicu     |   4 |                 3 | 697029           |     7438.27   |                             0.5521 |                    0.0311 |                             0.9999 |                   1.3788 |                               0 |                28.3371 |
| eicu     |   5 |                 3 | 698743           |     2967.86   |                             0.75   |                    0.0668 |                             0.9724 |                   1.3644 |                               0 |                34.3774 |
| mimic    |   2 |                 3 |      1.46479e+06 |       34.9525 |                             0      |                    0.0219 |                             0.9979 |                   0.0049 |                               1 |                38.3904 |
| mimic    |   3 |                 3 |      1.46315e+06 |      708.188  |                             0.0139 |                    0.0862 |                             0.968  |                   0      |                               2 |                53.6666 |
| mimic    |   4 |                 3 |      1.73394e+06 |   236084      |                             0.5    |                    0.0714 |                             0.9984 |                   1.3301 |                               0 |               103.798  |
| mimic    |   5 |                 3 |      1.88602e+06 |     1120.07   |                             0.7833 |                    0.0602 |                             0.9348 |                   1.4142 |                               0 |               116.825  |

Selection counts show how often each K minimized the archived two-axis score within a
cohort/seed grid. This fresh screen strengthens cross-cohort traceability but does not
retroactively validate the unsupported K=6–8 values in the submitted supplement.
The deeper archived eICU K=2–5 models remain the primary historical K-selection
evidence; disagreement in this screening grid must be reported, not overridden.
