#!/usr/bin/env python3
"""Summarize deep multi-seed eICU urine-output documentation sensitivities."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


SCENARIOS = ("documented_windows", "high_coverage")
SEEDS = (20260805, 20260806, 20260807)
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
        default=repo / "02_revision_outputs/intermediate/W3_uo_multiseed_refits",
    )
    parser.add_argument(
        "--original-data",
        type=Path,
        default=(
            repo
            / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
            / "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_uo_multiseed_sensitivity",
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


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    archived = pd.read_csv(args.original_data, usecols=["stay_id", "groupHPD"])
    archived = archived.drop_duplicates("stay_id")
    rows: list[dict[str, object]] = []
    missing: list[str] = []
    for scenario in SCENARIOS:
        for seed in SEEDS:
            stem = f"eicu__{scenario}_deep__K3__seed{seed}"
            assignment_path = args.refit_dir / f"{stem}__assignments.csv"
            diagnostic_path = args.refit_dir / f"{stem}__diagnostics.csv"
            if not assignment_path.exists() or not diagnostic_path.exists():
                missing.append(stem)
                continue
            assignments = pd.read_csv(assignment_path)
            diagnostics = pd.read_csv(diagnostic_path)
            merged = assignments.merge(archived, on="stay_id", validate="one_to_one")
            mapping = align_labels(merged.groupHPD, merged.group_median)
            aligned = merged.group_median.map(mapping)
            prevalence = aligned.value_counts(normalize=True).reindex([1, 2, 3], fill_value=0)
            agreement = float(aligned.eq(merged.groupHPD).mean())
            ari = float(adjusted_rand_score(merged.groupHPD, aligned))
            minimum_prevalence = float(prevalence.min())
            high_lag1 = float(diagnostics.loc[0, "high_absolute_lag1_fraction"])
            pass_fit = (
                agreement >= AGREEMENT_THRESHOLD
                and ari >= ARI_THRESHOLD
                and minimum_prevalence >= MINIMUM_CLUSTER_PREVALENCE
                and high_lag1 <= MAXIMUM_HIGH_LAG1_FRACTION
            )
            rows.append(
                {
                    "scenario": scenario,
                    "seed": seed,
                    "patients": len(merged),
                    "retained_fraction": len(merged) / len(archived),
                    "label_mapping": json.dumps(mapping, sort_keys=True),
                    "exact_agreement": agreement,
                    "ari": ari,
                    "nmi": float(
                        normalized_mutual_info_score(merged.groupHPD, aligned)
                    ),
                    "uncertain_fraction": float(
                        diagnostics.loc[0, "uncertain_fraction"]
                    ),
                    "minimum_cluster_prevalence": minimum_prevalence,
                    "dr_fraction": float(prevalence.loc[1]),
                    "rr_fraction": float(prevalence.loc[2]),
                    "pw_fraction": float(prevalence.loc[3]),
                    "high_absolute_lag1_fraction": high_lag1,
                    "fit_status": "PASS" if pass_fit else "CAUTION",
                }
            )
    if missing:
        status = {
            "overall_status": "FAIL_INCOMPLETE_UO_MULTI_SEED_GRID",
            "missing_refits": missing,
        }
        (args.output_dir / "uo_multiseed_status.json").write_text(
            json.dumps(status, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(status, indent=2))
        return 1

    metrics = pd.DataFrame(rows)
    summary = metrics.groupby("scenario", as_index=False).agg(
        seeds_completed=("seed", "nunique"),
        passing_seeds=("fit_status", lambda values: int(values.eq("PASS").sum())),
        patients=("patients", "first"),
        retained_fraction=("retained_fraction", "first"),
        median_exact_agreement=("exact_agreement", "median"),
        minimum_exact_agreement=("exact_agreement", "min"),
        median_ari=("ari", "median"),
        minimum_ari=("ari", "min"),
        median_nmi=("nmi", "median"),
        maximum_uncertain_fraction=("uncertain_fraction", "max"),
        minimum_cluster_prevalence=("minimum_cluster_prevalence", "min"),
        median_dr_fraction=("dr_fraction", "median"),
        median_rr_fraction=("rr_fraction", "median"),
        median_pw_fraction=("pw_fraction", "median"),
        maximum_high_lag1_fraction=("high_absolute_lag1_fraction", "max"),
    )
    summary["scenario_status"] = np.where(
        summary.passing_seeds.eq(3),
        "PASS_ALL_3_INITIALIZATIONS",
        "CAUTION_DOCUMENTATION_SELECTION_OR_MODE_SENSITIVITY",
    )
    status = {
        "overall_status": (
            "PASS_ALL_UO_MULTI_SEED_SCENARIOS"
            if summary.passing_seeds.eq(3).all()
            else "UO_MULTI_SEED_SENSITIVITY_CAUTIONS_REMAIN"
        ),
        "expected_fits": len(SCENARIOS) * len(SEEDS),
        "observed_fits": len(metrics),
        "all_prespecified_seeds_retained": True,
        "scenario_status": summary[["scenario", "scenario_status"]].to_dict(
            orient="records"
        ),
    }
    metrics.to_csv(args.output_dir / "uo_multiseed_seed_metrics.csv", index=False)
    summary.to_csv(args.output_dir / "uo_multiseed_scenario_summary.csv", index=False)
    (args.output_dir / "uo_multiseed_status.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    report = f"""# Deep multi-seed urine-output documentation sensitivity

**Status: {status['overall_status']}**

Both reviewer-requested eICU documentation scenarios were run with three deep
independent initializations. Every seed is retained. A scenario passes only when
all three fits satisfy the prespecified mixing, agreement, ARI, and minimum-cluster
prevalence thresholds.

{summary.round(4).to_markdown(index=False)}

## Seed-level evidence

{metrics.round(4).to_markdown(index=False)}

Restriction to documented or high-coverage windows changes the target population;
agreement in these selected subsets does not establish robustness for all 1,417
original eICU patients.
"""
    (args.output_dir / "W3_UO_MULTI_SEED_SENSITIVITY.md").write_text(
        report, encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
