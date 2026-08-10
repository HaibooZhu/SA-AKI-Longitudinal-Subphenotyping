#!/usr/bin/env python3
"""Summarize the complete fresh K=2-5, three-seed mixAK screening grid."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


COHORTS = ("mimic", "eicu", "aumc")
KS = (2, 3, 4, 5)
SEEDS = (20260805, 20260806, 20260807)


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=repo / "02_revision_outputs/intermediate/W3_fresh_refits",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_fresh_k_grid",
    )
    return parser.parse_args()


def expected_grid() -> set[tuple[str, int, int]]:
    return {(cohort, k, seed) for cohort in COHORTS for k in KS for seed in SEEDS}


def minmax(values: pd.Series) -> pd.Series:
    minimum, maximum = values.min(), values.max()
    if np.isclose(minimum, maximum):
        return pd.Series(0.0, index=values.index)
    return (values - minimum) / (maximum - minimum)


def score_grid(diagnostics: pd.DataFrame) -> pd.DataFrame:
    scored = diagnostics.copy()
    scored["scaled_deviance"] = scored.groupby(
        ["cohort", "seed"], group_keys=False
    )["mean_deviance"].transform(minmax)
    scored["scaled_lag1_failure"] = scored.groupby(
        ["cohort", "seed"], group_keys=False
    )["high_absolute_lag1_fraction"].transform(minmax)
    scored["historical_selection_score"] = np.sqrt(
        scored.scaled_deviance**2 + scored.scaled_lag1_failure**2
    )
    scored["rank_within_cohort_seed"] = scored.groupby(["cohort", "seed"])[
        "historical_selection_score"
    ].rank(method="min")
    scored["selected_within_seed"] = scored.rank_within_cohort_seed.eq(1)
    return scored


def summarize(scored: pd.DataFrame) -> pd.DataFrame:
    return (
        scored.groupby(["cohort", "K"], as_index=False)
        .agg(
            seeds_completed=("seed", "nunique"),
            mean_deviance=("mean_deviance", "mean"),
            sd_deviance=("mean_deviance", "std"),
            mean_high_absolute_lag1_fraction=(
                "high_absolute_lag1_fraction",
                "mean",
            ),
            mean_uncertain_fraction=("uncertain_fraction", "mean"),
            median_max_posterior_probability=(
                "median_max_posterior_probability",
                "median",
            ),
            median_selection_score=("historical_selection_score", "median"),
            selections_across_three_seeds=("selected_within_seed", "sum"),
            mean_elapsed_seconds=("elapsed_seconds", "mean"),
        )
        .sort_values(["cohort", "K"])
    )


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(args.input_dir.glob("*__primary__K*__diagnostics.csv"))
    diagnostics = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    observed = {
        (str(row.cohort), int(row.K), int(row.seed))
        for row in diagnostics.itertuples()
    }
    missing = sorted(expected_grid() - observed)
    duplicates = int(diagnostics.duplicated(["cohort", "K", "seed"]).sum())
    if missing or duplicates:
        status = {
            "overall_status": "FAIL_INCOMPLETE_FRESH_K_GRID",
            "expected_fits": len(expected_grid()),
            "observed_unique_fits": len(observed),
            "duplicate_fits": duplicates,
            "missing_fits": [f"{c}:K{k}:seed{s}" for c, k, s in missing],
        }
        (args.output_dir / "fresh_k_grid_status.json").write_text(
            json.dumps(status, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(status, indent=2))
        return 1

    scored = score_grid(diagnostics)
    summary = summarize(scored)
    selected = {
        cohort: {
            str(int(row.K)): int(row.selections_across_three_seeds)
            for row in summary.loc[summary.cohort.eq(cohort)].itertuples()
            if row.selections_across_three_seeds > 0
        }
        for cohort in COHORTS
    }
    unanimous_k3 = all(selected[cohort] == {"3": 3} for cohort in COHORTS)
    status = {
        "overall_status": (
            "PASS_COMPLETE_GRID_UNANIMOUS_K3"
            if unanimous_k3
            else "COMPLETE_GRID_NONUNANIMOUS_K3"
        ),
        "expected_fits": len(expected_grid()),
        "observed_unique_fits": len(observed),
        "candidate_k": list(KS),
        "seeds": list(SEEDS),
        "selection_counts": selected,
        "probability_scaling": "native 0-1",
        "screening_depth_limitation": (
            "Fresh models use the uniform screening-depth configuration. The deeper "
            "archived eICU K=2-5 audit remains the primary historical evidence."
        ),
    }
    scored.to_csv(args.output_dir / "fresh_k_grid_all_seed_diagnostics.csv", index=False)
    summary.to_csv(args.output_dir / "fresh_k_grid_summary.csv", index=False)
    (args.output_dir / "fresh_k_grid_status.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    report = f"""# Fresh cross-cohort K=2–5 screening grid

**Status: {status['overall_status']}**

All 36 prespecified fits were completed: three cohorts × four candidate K values ×
three deterministic seeds. Every fit used the same screening-depth burn, keep, thin,
lag-1 diagnostic, HPD uncertainty rule, and native 0–1 posterior probability scale.

{summary.round(4).to_markdown(index=False)}

Selection counts show how often each K minimized the archived two-axis score within a
cohort/seed grid. This fresh screen strengthens cross-cohort traceability but does not
retroactively validate the unsupported K=6–8 values in the submitted supplement.
The deeper archived eICU K=2–5 models remain the primary historical K-selection
evidence; disagreement in this screening grid must be reported, not overridden.
"""
    (args.output_dir / "W3_FRESH_CROSS_COHORT_K_GRID.md").write_text(
        report, encoding="utf-8"
    )
    print(summary.to_string(index=False))
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
