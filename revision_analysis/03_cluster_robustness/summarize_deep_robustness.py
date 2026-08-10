#!/usr/bin/env python3
"""Summarize deep three-seed reruns of all screening-depth CAUTION scenarios."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


SCENARIOS = (
    ("mimic", "complete_30_window_followup"),
    ("mimic", "limited_forward_fill_complete_rows"),
    ("eicu", "complete_30_window_followup"),
    ("eicu", "limited_forward_fill_complete_rows"),
    ("aumc", "exclude_documented_rrt"),
    ("aumc", "complete_30_window_followup"),
)
SEEDS = (20260805, 20260806, 20260807)
COHORT_PATHS = {
    "mimic": "01.MIMICIV_SAKI_trajCluster/df_mixAK_fea4_C3.csv",
    "eicu": "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv",
    "aumc": "02.AUMCdb_SAKI_trajCluster/df_mixAK_fea3_C3_aumc.csv",
}
AGREEMENT_THRESHOLD = 0.75
ARI_THRESHOLD = 0.50
MINIMUM_CLUSTER_PREVALENCE = 0.03
MAXIMUM_HIGH_LAG1_FRACTION = 0.0


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refit-dir",
        type=Path,
        default=repo / "02_revision_outputs/intermediate/W3_deep_robustness_refits",
    )
    parser.add_argument(
        "--snapshot-root",
        type=Path,
        default=repo / "00_frozen_inputs/data_snapshot/remote_project_snapshot",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_deep_robustness",
    )
    return parser.parse_args()


def align_labels(reference: pd.Series, candidate: pd.Series) -> dict[int, int]:
    contingency = pd.crosstab(candidate, reference).reindex(
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
    seed: int,
    assignments: pd.DataFrame,
    diagnostics: pd.DataFrame,
    archived: pd.DataFrame,
) -> dict[str, object]:
    merged = assignments.merge(
        archived[["stay_id", "groupHPD"]].drop_duplicates("stay_id"),
        on="stay_id",
        validate="one_to_one",
    )
    mapping = align_labels(merged.groupHPD, merged.group_median)
    aligned = merged.group_median.map(mapping)
    prevalence = aligned.value_counts(normalize=True).reindex([1, 2, 3], fill_value=0)
    agreement = float(aligned.eq(merged.groupHPD).mean())
    ari = float(adjusted_rand_score(merged.groupHPD, aligned))
    minimum_prevalence = float(prevalence.min())
    high_lag1 = float(diagnostics.loc[0, "high_absolute_lag1_fraction"])
    fit_pass = (
        agreement >= AGREEMENT_THRESHOLD
        and ari >= ARI_THRESHOLD
        and minimum_prevalence >= MINIMUM_CLUSTER_PREVALENCE
        and high_lag1 <= MAXIMUM_HIGH_LAG1_FRACTION
    )
    return {
        "cohort": cohort,
        "scenario": scenario,
        "seed": seed,
        "archived_patients": int(archived.stay_id.nunique()),
        "refit_patients": len(merged),
        "retained_fraction": len(merged) / archived.stay_id.nunique(),
        "label_mapping": json.dumps(mapping, sort_keys=True),
        "exact_agreement": agreement,
        "ari": ari,
        "nmi": float(normalized_mutual_info_score(merged.groupHPD, aligned)),
        "uncertain_fraction": float(diagnostics.loc[0, "uncertain_fraction"]),
        "minimum_cluster_prevalence": minimum_prevalence,
        "mean_deviance": float(diagnostics.loc[0, "mean_deviance"]),
        "high_absolute_lag1_fraction": high_lag1,
        "fit_status": "PASS" if fit_pass else "CAUTION",
    }


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    missing: list[str] = []
    archived_by_cohort = {
        cohort: pd.read_csv(args.snapshot_root / relative)
        for cohort, relative in COHORT_PATHS.items()
    }
    for cohort, scenario in SCENARIOS:
        for seed in SEEDS:
            stem = f"{cohort}__{scenario}_deep__K3__seed{seed}"
            assignment_path = args.refit_dir / f"{stem}__assignments.csv"
            diagnostic_path = args.refit_dir / f"{stem}__diagnostics.csv"
            if not assignment_path.exists() or not diagnostic_path.exists():
                missing.append(stem)
                continue
            diagnostics = pd.read_csv(diagnostic_path)
            if len(diagnostics) != 1:
                raise ValueError(f"{stem}: expected exactly one diagnostic row")
            rows.append(
                evaluate(
                    cohort,
                    scenario,
                    seed,
                    pd.read_csv(assignment_path),
                    diagnostics,
                    archived_by_cohort[cohort],
                )
            )
    if missing:
        status = {
            "overall_status": "FAIL_INCOMPLETE_DEEP_CAUTION_GRID",
            "missing_refits": missing,
        }
        (args.output_dir / "deep_robustness_status.json").write_text(
            json.dumps(status, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(status, indent=2))
        return 1

    metrics = pd.DataFrame(rows)
    summary = (
        metrics.groupby(["cohort", "scenario"], as_index=False)
        .agg(
            seeds_completed=("seed", "nunique"),
            passing_seeds=("fit_status", lambda values: int(values.eq("PASS").sum())),
            median_exact_agreement=("exact_agreement", "median"),
            minimum_exact_agreement=("exact_agreement", "min"),
            median_ari=("ari", "median"),
            minimum_ari=("ari", "min"),
            median_nmi=("nmi", "median"),
            maximum_uncertain_fraction=("uncertain_fraction", "max"),
            minimum_cluster_prevalence=("minimum_cluster_prevalence", "min"),
            maximum_high_lag1_fraction=("high_absolute_lag1_fraction", "max"),
        )
    )
    summary["scenario_status"] = np.where(
        summary.passing_seeds.eq(3),
        "PASS_ALL_3_INITIALIZATIONS",
        "CAUTION_MODE_OR_DATA_SENSITIVITY",
    )
    cautions = summary.loc[
        summary.scenario_status.ne("PASS_ALL_3_INITIALIZATIONS"),
        ["cohort", "scenario", "scenario_status"],
    ]
    status = {
        "overall_status": (
            "PASS_ALL_DEEP_CAUTION_RERUNS"
            if cautions.empty
            else "DEEP_RERUNS_RETAIN_SENSITIVITY_CAUTIONS"
        ),
        "expected_fits": len(SCENARIOS) * len(SEEDS),
        "observed_fits": len(metrics),
        "thresholds": {
            "exact_agreement_minimum": AGREEMENT_THRESHOLD,
            "ari_minimum": ARI_THRESHOLD,
            "minimum_cluster_prevalence": MINIMUM_CLUSTER_PREVALENCE,
            "maximum_high_lag1_fraction": MAXIMUM_HIGH_LAG1_FRACTION,
        },
        "caution_scenarios": cautions.to_dict(orient="records"),
        "all_prespecified_seeds_retained": True,
    }
    metrics.to_csv(args.output_dir / "deep_robustness_seed_metrics.csv", index=False)
    summary.to_csv(args.output_dir / "deep_robustness_scenario_summary.csv", index=False)
    (args.output_dir / "deep_robustness_status.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    report = f"""# Deep reruns of screening-depth CAUTION scenarios

**Status: {status['overall_status']}**

All 18 prespecified fits were retained: six cautionary cohort/scenario pairs ×
three deep independent initializations. A scenario passes only when all three
initializations meet the agreement, ARI, minimum-prevalence, and lag-1 thresholds.
Failed or degenerate fits are not removed from the denominator.

{summary.round(4).to_markdown(index=False)}

## Seed-level evidence

{metrics.round(4).to_markdown(index=False)}

These reruns distinguish persistent data-processing sensitivity from a favorable
single short-chain result; they do not convert a caution to pass by visual inspection.
"""
    (args.output_dir / "W3_DEEP_CAUTION_RERUNS.md").write_text(
        report, encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
