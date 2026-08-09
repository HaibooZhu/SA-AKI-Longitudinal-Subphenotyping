# Legacy baseline limitations and revision authority

## Status

The `src/`, `scripts/`, and `configs/` directories preserve an earlier public
development snapshot. They are retained unchanged for provenance and are not the
authoritative executable basis for the JTIM revision. Corrected analyses are isolated
under `revision_analysis/`.

## Confirmed legacy limitations

| Legacy location | Confirmed issue | Revision replacement |
| --- | --- | --- |
| `src/sa_aki_pipeline/preprocessing/time_windows.py` | The historical grouped forward-fill path can drop the grouping identifier from the returned frame. | `revision_analysis/02_missingness_sensitivity/build_corrected_time_grid.py` preserves identifiers and asserts a complete patient × 30-window grid. |
| `src/sa_aki_pipeline/phenotyping/mixak.py` | The historical export divides posterior medians by two for display and constructs the random-intercept vector from a scalar K rather than the outcome count. | Revision mixAK runners retain native 0–1 probabilities and set one random-intercept flag per renal feature. |
| `src/sa_aki_pipeline/fluid/config.py` and historical YAML templates | Historical diuretic/urine-output settings are not a verified methods contract and conflict with the traced revision definitions. | Revision workstreams explicitly audit the historical treatment analysis and reconstruct the eICU six-hour urine-output aggregation from frozen raw events. |

## Reporting rule

No numerical result, methods statement, or reviewer response should cite successful
execution of the legacy baseline package as validation. Revision claims must point to
the corresponding `revision_analysis/` script, its aggregate audit output, and the
exact Git commit. This separation keeps the original code fixed while preventing it
from being mistaken for corrected software.
