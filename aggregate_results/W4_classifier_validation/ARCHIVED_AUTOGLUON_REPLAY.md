# Archived AutoGluon classifier replay

## Status

The submitted `XGBoost_BAG_L2` model was successfully reloaded and evaluated using the
archived software versions: AutoGluon 0.7.0, scikit-learn
1.2.2, and XGBoost 1.7.4. Only aggregate outputs
were exported; no patient identifiers or patient-level predictions were written.

## Primary model and simple comparator

| model                             | cohort              |    n |   accuracy |   balanced_accuracy |   precision_macro |   recall_macro |   f1_macro |   auc_ovo_macro |   auc_ovr_macro |   log_loss |   multiclass_brier |
|:----------------------------------|:--------------------|-----:|-----------:|--------------------:|------------------:|---------------:|-----------:|----------------:|----------------:|-----------:|-------------------:|
| archived_AutoGluon_XGBoost_BAG_L2 | internal_MIMIC_eICU | 1227 |      0.792 |               0.637 |             0.745 |          0.637 |      0.674 |           0.836 |           0.869 |      0.536 |              0.308 |
| simple_six_variable_logistic      | internal_MIMIC_eICU | 1227 |      0.768 |               0.6   |             0.721 |          0.6   |      0.639 |           0.801 |           0.836 |      0.579 |              0.337 |
| archived_AutoGluon_XGBoost_BAG_L2 | external_AUMC       | 2183 |      0.555 |               0.528 |             0.56  |          0.528 |      0.504 |           0.769 |           0.789 |      0.92  |              0.576 |
| simple_six_variable_logistic      | external_AUMC       | 2183 |      0.713 |               0.578 |             0.599 |          0.578 |      0.583 |           0.769 |           0.805 |      0.697 |              0.407 |

The comparator is a prespecified multinomial logistic regression using six transparent
first-24-hour kidney summaries: creatinine minimum/maximum, creatinine-to-baseline ratio
minimum/maximum, and urine-output minimum/mean. It is a deliberately simple reference,
not a proposed clinical score.

## One-versus-rest calibration

| model                             | cohort              |   class | phenotype   |   calibration_intercept |   calibration_slope |   ece_deciles |
|:----------------------------------|:--------------------|--------:|:------------|------------------------:|--------------------:|--------------:|
| archived_AutoGluon_XGBoost_BAG_L2 | internal_MIMIC_eICU |       1 | DR          |                  -0.098 |               0.821 |         0.04  |
| archived_AutoGluon_XGBoost_BAG_L2 | internal_MIMIC_eICU |       2 | RR          |                   0.036 |               0.976 |         0.024 |
| archived_AutoGluon_XGBoost_BAG_L2 | internal_MIMIC_eICU |       3 | PW          |                   0.047 |               1.092 |         0.017 |
| simple_six_variable_logistic      | internal_MIMIC_eICU |       1 | DR          |                   0.054 |               1.061 |         0.04  |
| simple_six_variable_logistic      | internal_MIMIC_eICU |       2 | RR          |                  -0.017 |               1.034 |         0.041 |
| simple_six_variable_logistic      | internal_MIMIC_eICU |       3 | PW          |                   0.021 |               1.04  |         0.016 |
| archived_AutoGluon_XGBoost_BAG_L2 | external_AUMC       |       1 | DR          |                  -1.1   |               0.437 |         0.261 |
| archived_AutoGluon_XGBoost_BAG_L2 | external_AUMC       |       2 | RR          |                   1.866 |               1.004 |         0.277 |
| archived_AutoGluon_XGBoost_BAG_L2 | external_AUMC       |       3 | PW          |                  -0.831 |               0.784 |         0.049 |
| simple_six_variable_logistic      | external_AUMC       |       1 | DR          |                  -0.584 |               0.855 |         0.113 |
| simple_six_variable_logistic      | external_AUMC       |       2 | RR          |                   0.788 |               0.931 |         0.126 |
| simple_six_variable_logistic      | external_AUMC       |       3 | PW          |                  -0.766 |               0.683 |         0.025 |

Ideal calibration has intercept 0 and slope 1. ECE is an absolute decile-weighted
calibration error; lower is better.

## Pair-normalized discrimination

| model                             | cohort              | comparison   |   n_pair |   auc_pair_normalized |
|:----------------------------------|:--------------------|:-------------|---------:|----------------------:|
| archived_AutoGluon_XGBoost_BAG_L2 | internal_MIMIC_eICU | DR vs RR     |     1137 |                 0.863 |
| archived_AutoGluon_XGBoost_BAG_L2 | internal_MIMIC_eICU | DR vs PW     |      429 |                 0.774 |
| archived_AutoGluon_XGBoost_BAG_L2 | internal_MIMIC_eICU | RR vs PW     |      888 |                 0.947 |
| simple_six_variable_logistic      | internal_MIMIC_eICU | DR vs RR     |     1137 |                 0.817 |
| simple_six_variable_logistic      | internal_MIMIC_eICU | DR vs PW     |      429 |                 0.791 |
| simple_six_variable_logistic      | internal_MIMIC_eICU | RR vs PW     |      888 |                 0.929 |
| archived_AutoGluon_XGBoost_BAG_L2 | external_AUMC       | DR vs RR     |     1964 |                 0.818 |
| archived_AutoGluon_XGBoost_BAG_L2 | external_AUMC       | DR vs PW     |      783 |                 0.653 |
| archived_AutoGluon_XGBoost_BAG_L2 | external_AUMC       | RR vs PW     |     1619 |                 0.933 |
| simple_six_variable_logistic      | external_AUMC       | DR vs RR     |     1964 |                 0.828 |
| simple_six_variable_logistic      | external_AUMC       | DR vs PW     |      783 |                 0.665 |
| simple_six_variable_logistic      | external_AUMC       | RR vs PW     |     1619 |                 0.913 |

These are the three undirected phenotype pairs. For each pair, the score is
`p(class A) / [p(class A) + p(class B)]`; therefore no six-direction raw-probability
results are mislabeled as ordinary one-versus-one AUCs.

## Paired incremental-value analysis

| cohort              | metric            |   archived_autogluon |   simple_logistic |   difference_autogluon_minus_simple |   ci_low |   ci_high | direction        |   bootstrap_replicates |
|:--------------------|:------------------|---------------------:|------------------:|------------------------------------:|---------:|----------:|:-----------------|-----------------------:|
| internal_MIMIC_eICU | balanced_accuracy |                0.637 |             0.6   |                               0.037 |    0.007 |     0.065 | higher_is_better |                   1000 |
| internal_MIMIC_eICU | f1_macro          |                0.674 |             0.639 |                               0.035 |    0.003 |     0.066 | higher_is_better |                   1000 |
| internal_MIMIC_eICU | auc_ovo_macro     |                0.836 |             0.801 |                               0.035 |    0.022 |     0.048 | higher_is_better |                   1000 |
| internal_MIMIC_eICU | multiclass_brier  |                0.308 |             0.337 |                              -0.029 |   -0.042 |    -0.015 | lower_is_better  |                   1000 |
| external_AUMC       | balanced_accuracy |                0.528 |             0.578 |                              -0.05  |   -0.071 |    -0.028 | higher_is_better |                   1000 |
| external_AUMC       | f1_macro          |                0.504 |             0.583 |                              -0.079 |   -0.101 |    -0.058 | higher_is_better |                   1000 |
| external_AUMC       | auc_ovo_macro     |                0.769 |             0.769 |                              -0     |   -0.013 |     0.012 | higher_is_better |                   1000 |
| external_AUMC       | multiclass_brier  |                0.576 |             0.407 |                               0.169 |    0.155 |     0.183 | lower_is_better  |                   1000 |

Differences are AutoGluon minus simple logistic on the same test patients with
class-stratified bootstrap intervals. Positive differences favor AutoGluon for AUC,
F1, and balanced accuracy; negative differences favor AutoGluon for Brier score.
