#!/usr/bin/env python3
"""Align K=3 robustness refits to archived labels and export aggregate metrics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


SEED = 20260805
SCENARIOS = (
    "zscore_scaling",
    "robust_scaling",
    "exclude_documented_rrt",
    "complete_30_window_followup",
    "limited_forward_fill_complete_rows",
)
COHORTS = {
    "mimic": "01.MIMICIV_SAKI_trajCluster/df_mixAK_fea4_C3.csv",
    "eicu": "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv",
    "aumc": "02.AUMCdb_SAKI_trajCluster/df_mixAK_fea3_C3_aumc.csv",
}
AGREEMENT_THRESHOLD = 0.75
ARI_THRESHOLD = 0.50
MINIMUM_CLUSTER_PREVALENCE = 0.03


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refit-dir",
        type=Path,
        default=repo / "02_revision_outputs/intermediate/W3_robustness_refits",
    )
    parser.add_argument(
        "--snapshot-root",
        type=Path,
        default=repo / "00_frozen_inputs/data_snapshot/remote_project_snapshot",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_cross_cohort_robustness",
    )
    return parser.parse_args()


def align_labels(original: pd.Series, refit: pd.Series) -> dict[int, int]:
    contingency = pd.crosstab(refit, original).reindex(
        index=[1, 2, 3], columns=[1, 2, 3], fill_value=0
    )
    rows, columns = linear_sum_assignment(-contingency.to_numpy())
    return {
        int(contingency.index[row]): int(contingency.columns[column])
        for row, column in zip(rows, columns)
    }


def evaluate(
    cohort: str,
    scenario: str,
    assignments: pd.DataFrame,
    archived: pd.DataFrame,
    diagnostics: pd.DataFrame,
) -> dict[str, object]:
    labels = archived[["stay_id", "groupHPD"]].drop_duplicates("stay_id")
    merged = assignments.merge(labels, on="stay_id", how="inner", validate="one_to_one")
    mapping = align_labels(merged.groupHPD, merged.group_median)
    merged["aligned_group"] = merged.group_median.map(mapping)
    merged["aligned_hpd"] = merged.group_hpd.map(mapping).fillna(4).astype(int)
    confident = merged.loc[merged.aligned_hpd.ne(4)]
    prevalence = merged.aligned_group.value_counts(normalize=True).reindex([1, 2, 3], fill_value=0)
    agreement = float(merged.aligned_group.eq(merged.groupHPD).mean())
    ari = float(adjusted_rand_score(merged.groupHPD, merged.aligned_group))
    minimum_prevalence = float(prevalence.min())
    return {
        "cohort": cohort,
        "scenario": scenario,
        "seed": SEED,
        "archived_patients": int(labels.stay_id.nunique()),
        "refit_patients": int(len(merged)),
        "retained_fraction": float(len(merged) / labels.stay_id.nunique()),
        "label_mapping": json.dumps(mapping, sort_keys=True),
        "exact_agreement": agreement,
        "ari": ari,
        "nmi": float(normalized_mutual_info_score(merged.groupHPD, merged.aligned_group)),
        "uncertain_fraction": float(merged.aligned_hpd.eq(4).mean()),
        "confident_patients": int(len(confident)),
        "confident_agreement": (
            float(confident.aligned_hpd.eq(confident.groupHPD).mean())
            if len(confident)
            else np.nan
        ),
        "minimum_cluster_prevalence": minimum_prevalence,
        "mean_deviance": float(diagnostics.loc[0, "mean_deviance"]),
        "high_absolute_lag1_fraction": float(
            diagnostics.loc[0, "high_absolute_lag1_fraction"]
        ),
        "screening_threshold_status": (
            "PASS"
            if agreement >= AGREEMENT_THRESHOLD
            and ari >= ARI_THRESHOLD
            and minimum_prevalence >= MINIMUM_CLUSTER_PREVALENCE
            else "CAUTION"
        ),
    }


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    missing: list[str] = []
    for cohort, relative_matrix in COHORTS.items():
        archived = pd.read_csv(args.snapshot_root / relative_matrix)
        for scenario in SCENARIOS:
            stem = f"{cohort}__{scenario}__K3__seed{SEED}"
            assignment_path = args.refit_dir / f"{stem}__assignments.csv"
            diagnostic_path = args.refit_dir / f"{stem}__diagnostics.csv"
            if not assignment_path.exists() or not diagnostic_path.exists():
                missing.append(stem)
                continue
            diagnostics = pd.read_csv(diagnostic_path)
            if len(diagnostics) != 1:
                raise ValueError(f"{stem}: expected one diagnostic row")
            rows.append(
                evaluate(
                    cohort,
                    scenario,
                    pd.read_csv(assignment_path),
                    archived,
                    diagnostics,
                )
            )
    if missing:
        status = {
            "overall_status": "FAIL_INCOMPLETE_ROBUSTNESS_REFITS",
            "missing_refits": missing,
        }
        (args.output_dir / "robustness_refit_status.json").write_text(
            json.dumps(status, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(status, indent=2))
        return 1

    metrics = pd.DataFrame(rows)
    caution = metrics.loc[
        metrics.screening_threshold_status.ne("PASS"), ["cohort", "scenario"]
    ]
    status = {
        "overall_status": (
            "PASS_ALL_PRESPECIFIED_ROBUSTNESS_THRESHOLDS"
            if caution.empty
            else "COMPLETE_WITH_SENSITIVITY_CAUTIONS"
        ),
        "completed_refits": len(metrics),
        "thresholds": {
            "exact_agreement_minimum": AGREEMENT_THRESHOLD,
            "ari_minimum": ARI_THRESHOLD,
            "minimum_cluster_prevalence": MINIMUM_CLUSTER_PREVALENCE,
        },
        "caution_scenarios": caution.to_dict(orient="records"),
        "post_rrt_truncation_available_across_all_cohorts": False,
    }
    metrics.to_csv(args.output_dir / "cross_cohort_robustness_metrics.csv", index=False)
    (args.output_dir / "robustness_refit_status.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    report = f"""# Cross-cohort K=3 robustness refits

**Status: {status['overall_status']}**

Numeric mixture labels were aligned to the archived phenotype labels by maximum
overlap before agreement, ARI, and NMI were calculated. Thresholds were prespecified
as exact agreement ≥ {AGREEMENT_THRESHOLD:.2f}, ARI ≥ {ARI_THRESHOLD:.2f}, and every
aligned cluster prevalence ≥ {MINIMUM_CLUSTER_PREVALENCE:.2f}.

{metrics.round(4).to_markdown(index=False)}

`CAUTION` is retained in the report and rebuttal; it is not converted to PASS based on
the visual similarity of trajectories. The complete-follow-up subset is selection
sensitive. The RRT analysis excludes documented recipients, while consistent
post-initiation truncation is unavailable because exact RRT start time is not present
for all three frozen cohorts. The limited-forward-fill scenario retains only renal
feature-complete rows and is the prespecified alternative to multiple imputation.
"""
    (args.output_dir / "W3_CROSS_COHORT_ROBUSTNESS_REFITS.md").write_text(
        report, encoding="utf-8"
    )
    print(metrics.to_string(index=False))
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
