# W6 Diuretic analysis audit and exploratory landmark sensitivity analysis

## Bottom line

The archived diuretic analysis cannot identify a treatment effect. It classifies
patients using urine output after treatment, combines responses across an unrestricted
number of post-onset administrations using an ad hoc majority/tie rule, and compares
those post-treatment groups with patients who did not receive a diuretic. The archived
three-group propensity matching used only day-1 creatinine, urine output, nonrenal SOFA,
and colloid input and did not retain matched-triplet identifiers. The diuretic findings
must therefore be demoted to exploratory, associative supplementary results.

## Historical matching implementation audit

| cohort   | phenotype   | historical_method   | historical_tuning   | historical_covariates                                | source_file                        |
|:---------|:------------|:--------------------|:--------------------|:-----------------------------------------------------|:-----------------------------------|
| MIMIC-IV | DR          | TriMatch            | caliper=0.05        | creatinine, urineoutput, sofa_norenal, colloid_bolus | 04.R_diuretic_responsitive-3psm.py |
| MIMIC-IV | RR          | TriMatch OneToN     | M1=1.5; M2=4        | creatinine, urineoutput, sofa_norenal, colloid_bolus | 04.R_diuretic_responsitive-3psm.py |
| MIMIC-IV | PW          | TriMatch            | caliper=0.14        | creatinine, urineoutput, sofa_norenal, colloid_bolus | 04.R_diuretic_responsitive-3psm.py |
| AUMC     | DR          | TriMatch            | caliper=0.10        | creatinine, urineoutput, sofa_norenal, colloid_bolus | 04.R_diuretic_responsitive-3psm.py |
| AUMC     | RR          | TriMatch            | caliper=0.03        | creatinine, urineoutput, sofa_norenal, colloid_bolus | 04.R_diuretic_responsitive-3psm.py |
| AUMC     | PW          | TriMatch            | caliper=0.20        | creatinine, urineoutput, sofa_norenal, colloid_bolus | 04.R_diuretic_responsitive-3psm.py |

The frozen executable analysis used `TriMatch`, with phenotype- and cohort-specific
calipers or OneToN settings. This does not match the submitted Methods description of
a single binary logistic nearest-neighbor procedure with a 0.2-logit caliper. The
historical matched analysis is therefore retained only as an audited legacy result; it
is not used to estimate a revised causal treatment effect. Source hashes are available
in `historical_matching_spec.csv`.

## Archived matching balance audit

Maximum pairwise absolute standardized mean differences (SMDs) were calculated among
the three archived groups within each phenotype, before and after archived matching.
An SMD above 0.10 indicates residual imbalance that should not be dismissed by a
nonsignificant p value.

| cohort_label   | stage           | used_in_archived_matching   |   max_smd |   median_smd |   variables |
|:---------------|:----------------|:----------------------------|----------:|-------------:|------------:|
| AUMC           | after_matching  | False                       |     0.551 |        0.278 |          12 |
| AUMC           | after_matching  | True                        |     0.298 |        0.138 |          12 |
| AUMC           | before_matching | False                       |     0.542 |        0.264 |          12 |
| AUMC           | before_matching | True                        |     0.653 |        0.315 |          12 |
| MIMIC-IV       | after_matching  | False                       |     0.432 |        0.187 |          12 |
| MIMIC-IV       | after_matching  | True                        |     0.387 |        0.151 |          12 |
| MIMIC-IV       | before_matching | False                       |     0.455 |        0.25  |          12 |
| MIMIC-IV       | before_matching | True                        |     0.613 |        0.195 |          12 |

Full variable-level results are in `archived_psm_balance_smd.csv`. Variables not used
in the archived matching include age, sex, weight, and baseline creatinine. Other
important unavailable or unaligned confounders include treatment indication, clinician
assessment of congestion, pre-onset fluid exposure, contemporaneous vasopressor dose,
nephrotoxin exposure, dose route, treatment limitation, and the exact timing of all
covariates relative to the first dose.

## Restricted 24-hour landmark cohort

To reduce, but not eliminate, immortal-time and exposure-timing bias, the sensitivity
analysis included only patients whose first post-SA-AKI diuretic administration occurred
in day 1, who were alive beyond 24 hours, and who had positive recorded urine volumes in
both two-hour windows around the first administration. The comparison is response versus
nonresponse among treated patients; untreated patients are not used as a causal control.

The archived first-dose rule labels a dose responsive if post-dose two-hour urine output
is greater than 200 mL when pre-dose output is already greater than 200 mL, or if output
is both greater than 200 mL and more than 10% above the pre-dose value. This makes the
archived rule nearly equivalent to the absolute 200-mL rule in these data. The stricter
10%-increase-plus-200-mL rule is reported separately.

| cohort_label   | phenotype   | response      |   n |   deaths |   mortality_risk |
|:---------------|:------------|:--------------|----:|---------:|-----------------:|
| MIMIC-IV       | DR          | nonresponsive | 107 |       16 |            0.15  |
| MIMIC-IV       | DR          | responsive    | 168 |       24 |            0.143 |
| MIMIC-IV       | RR          | nonresponsive | 126 |        9 |            0.071 |
| MIMIC-IV       | RR          | responsive    | 651 |       12 |            0.018 |
| MIMIC-IV       | PW          | nonresponsive |  48 |       21 |            0.438 |
| MIMIC-IV       | PW          | responsive    |  24 |       10 |            0.417 |
| AUMC           | DR          | nonresponsive |  82 |       23 |            0.28  |
| AUMC           | DR          | responsive    |  48 |       10 |            0.208 |
| AUMC           | RR          | nonresponsive |  56 |        8 |            0.143 |
| AUMC           | RR          | responsive    | 172 |       24 |            0.14  |
| AUMC           | PW          | nonresponsive |  22 |        7 |            0.318 |
| AUMC           | PW          | responsive    |  12 |        6 |            0.5   |

## Pooled adjusted response association under alternate response rules

Models adjust for phenotype, age, sex, baseline creatinine, day-1 nonrenal SOFA,
pre-dose two-hour urine output, and first-dose amount. The original first-dose rule is
shown alongside stricter alternatives to expose threshold dependence.

| cohort_label   | response_definition       |   adjusted_or_response_vs_nonresponse |   ci_low |   ci_high |   p_value |   n_complete |   deaths |   responsive_n |   n_parameters |   events_per_parameter | converged   | diagnostic_flag    |
|:---------------|:--------------------------|--------------------------------------:|---------:|----------:|----------:|-------------:|---------:|---------------:|---------------:|-----------------------:|:------------|:-------------------|
| MIMIC-IV       | response_archived_first   |                                 0.592 |    0.331 |     1.06  |     0.078 |         1124 |       92 |            843 |             10 |                    9.2 | True        | CAUTION_LOW_EVENTS |
| MIMIC-IV       | response_strict_10pct_200 |                                 0.601 |    0.359 |     1.006 |     0.053 |         1124 |       92 |            801 |             10 |                    9.2 | True        | CAUTION_LOW_EVENTS |
| MIMIC-IV       | response_absolute_200     |                                 0.592 |    0.331 |     1.06  |     0.078 |         1124 |       92 |            843 |             10 |                    9.2 | True        | CAUTION_LOW_EVENTS |
| AUMC           | response_archived_first   |                                 1.129 |    0.567 |     2.249 |     0.729 |          382 |       76 |            226 |             10 |                    7.6 | True        | CAUTION_LOW_EVENTS |
| AUMC           | response_strict_10pct_200 |                                 1.205 |    0.651 |     2.23  |     0.552 |          382 |       76 |            179 |             10 |                    7.6 | True        | CAUTION_LOW_EVENTS |
| AUMC           | response_absolute_200     |                                 1.117 |    0.559 |     2.229 |     0.755 |          382 |       76 |            227 |             10 |                    7.6 | True        | CAUTION_LOW_EVENTS |

## Phenotype-specific exploratory associations

| cohort_label   | phenotype   |   adjusted_or_response_vs_nonresponse |   ci_low |   ci_high |   p_value |   n_complete |   deaths |   n_parameters |   events_per_parameter | converged   | diagnostic_flag    |
|:---------------|:------------|--------------------------------------:|---------:|----------:|----------:|-------------:|---------:|---------------:|-----------------------:|:------------|:-------------------|
| MIMIC-IV       | DR          |                                 0.959 |    0.456 |     2.017 |     0.913 |         1124 |       92 |             12 |                  7.667 | True        | CAUTION_LOW_EVENTS |
| MIMIC-IV       | RR          |                                 0.225 |    0.093 |     0.545 |     0.001 |         1124 |       92 |             12 |                  7.667 | True        | CAUTION_LOW_EVENTS |
| MIMIC-IV       | PW          |                                 0.754 |    0.248 |     2.294 |     0.619 |         1124 |       92 |             12 |                  7.667 | True        | CAUTION_LOW_EVENTS |
| AUMC           | DR          |                                 1.137 |    0.374 |     3.458 |     0.82  |          382 |       76 |             12 |                  6.333 | True        | CAUTION_LOW_EVENTS |
| AUMC           | RR          |                                 0.995 |    0.415 |     2.384 |     0.99  |          382 |       76 |             12 |                  6.333 | True        | CAUTION_LOW_EVENTS |
| AUMC           | PW          |                                 1.871 |    0.375 |     9.335 |     0.445 |          382 |       76 |             12 |                  6.333 | True        | CAUTION_LOW_EVENTS |

Joint interaction tests:

| cohort_label   | test                                     |   p_value |
|:---------------|:-----------------------------------------|----------:|
| MIMIC-IV       | response-by-phenotype interaction (2 df) |     0.04  |
| AUMC           | response-by-phenotype interaction (2 df) |     0.785 |

The nominal MIMIC-IV interaction was not reproduced in AUMC. Because phenotype is
defined using the same post-onset trajectory during which treatment response occurs,
these estimates are subject to post-exposure conditioning and cannot be interpreted as
treatment-effect modification. `CAUTION_LOW_EVENTS` denotes fewer than 10 deaths per
fitted coefficient.

## Interpretation and reporting decision

1. Remove diuretic responsiveness from the title, abstract conclusion, and main claim.
2. Do not use the no-diuretic group to estimate benefit or harm.
3. Present the restricted first-dose landmark analysis only in the supplement, labeled
   post-treatment prognostic association.
4. State that a standardized prospective furosemide stress protocol, pre-treatment
   eligibility criteria, treatment indication, and time-varying confounding control
   would be required to evaluate treatment effect or effect modification.
5. Do not recommend phenotype-guided diuretic treatment from these data.
