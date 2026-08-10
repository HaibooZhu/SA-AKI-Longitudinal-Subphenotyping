#!/usr/bin/env python3
"""Build eICU clustering inputs that do not treat undocumented urine output as zero."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


# The archived trajectory definition spans two pre-onset windows and 28
# post-onset windows. There is intentionally no time-zero label because the
# historical binning maps [0, 6) hours to time 1.
PLANNED_WINDOWS = (-2, -1, *range(1, 29))
EXACT_MATCH_THRESHOLD = 0.999
MAE_THRESHOLD = 0.01
DOMINANCE_MARGIN = 0.10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", type=Path, required=True)
    parser.add_argument("--onsets", type=Path, required=True)
    parser.add_argument("--cluster-input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def collect_documented_windows(events_path: Path, included_ids: set[int]) -> pd.DataFrame:
    frames = []
    for chunk in pd.read_csv(
        events_path,
        usecols=["stay_id", "charttime", "urineoutput"],
        chunksize=500_000,
    ):
        chunk = chunk[chunk["stay_id"].isin(included_ids) & chunk["urineoutput"].notna()]
        if not chunk.empty:
            frames.append(chunk)
    if not frames:
        raise RuntimeError("No documented urine-output records were found for included patients.")
    return pd.concat(frames, ignore_index=True)


def calculate_coverage(audited: pd.DataFrame) -> pd.DataFrame:
    """Calculate both available-row and fixed planned-window coverage.

    The high-coverage sensitivity uses the fixed 30-window denominator. The
    available-row fraction is retained only to document why the prior 743-patient
    selection was inflated when follow-up rows ended early.
    """
    coverage = (
        audited.groupby(["stay_id", "groupHPD"], as_index=False)
        .agg(
            available_cluster_windows=("time", "size"),
            documented_windows=("uo_documented", "sum"),
            documented_zero_windows=("documented_zero", "sum"),
        )
    )
    coverage["planned_windows"] = len(PLANNED_WINDOWS)
    coverage["documented_fraction_available"] = (
        coverage["documented_windows"] / coverage["available_cluster_windows"]
    )
    coverage["documented_fraction_planned"] = (
        coverage["documented_windows"] / coverage["planned_windows"]
    )
    return coverage


def aggregation_metrics(archived: pd.Series, reconstructed: pd.Series) -> dict[str, float]:
    """Return aggregate-only reconstruction diagnostics for one candidate rule."""
    archived_numeric = pd.to_numeric(archived, errors="coerce")
    reconstructed_numeric = pd.to_numeric(reconstructed, errors="coerce")
    evaluable = archived_numeric.notna() & reconstructed_numeric.notna()
    if not evaluable.any():
        return {
            "evaluable_windows": 0,
            "exact_match_fraction": 0.0,
            "mae": float("inf"),
            "median_absolute_error": float("inf"),
            "maximum_absolute_error": float("inf"),
        }
    error = (
        archived_numeric.loc[evaluable] - reconstructed_numeric.loc[evaluable]
    ).abs()
    exact = np.isclose(
        archived_numeric.loc[evaluable],
        reconstructed_numeric.loc[evaluable],
        rtol=1e-9,
        atol=1e-9,
    )
    return {
        "evaluable_windows": int(evaluable.sum()),
        "exact_match_fraction": float(exact.mean()),
        "mae": float(error.mean()),
        "median_absolute_error": float(error.median()),
        "maximum_absolute_error": float(error.max()),
    }


def audit_aggregation_rule(
    archived: pd.Series,
    reconstructed_sum: pd.Series,
    reconstructed_mean: pd.Series,
) -> dict[str, object]:
    """Select sum or mean only when the archived transformation is reconstructed.

    The decision is deliberately fail-closed. A candidate must reproduce at least
    99.9% of documented windows, have essentially zero median error and MAE no
    greater than 0.01 mL, and materially outperform the alternative rule.
    """
    candidates = {
        "sum": aggregation_metrics(archived, reconstructed_sum),
        "mean": aggregation_metrics(archived, reconstructed_mean),
    }
    ranked = sorted(
        candidates,
        key=lambda name: (
            -candidates[name]["exact_match_fraction"],
            candidates[name]["mae"],
        ),
    )
    winner, runner_up = ranked
    winning = candidates[winner]
    dominance = (
        winning["exact_match_fraction"]
        - candidates[runner_up]["exact_match_fraction"]
    )
    passed = (
        winning["exact_match_fraction"] >= EXACT_MATCH_THRESHOLD
        and winning["median_absolute_error"] <= 1e-9
        and winning["mae"] <= MAE_THRESHOLD
        and dominance >= DOMINANCE_MARGIN
    )
    discrepancy_windows = int(
        round(
            winning["evaluable_windows"]
            * (1.0 - winning["exact_match_fraction"])
        )
    )
    return {
        "status": (
            f"PASS_{winner.upper()}_RECONSTRUCTED"
            + (
                "_WITH_ISOLATED_SOURCE_DISCREPANCY"
                if discrepancy_windows > 0
                else ""
            )
            if passed
            else "FAIL_UO_AGGREGATION_NOT_RECONSTRUCTED"
        ),
        "passed": passed,
        "selected_rule": winner if passed else None,
        "discrepancy_windows": discrepancy_windows,
        "exact_match_threshold": EXACT_MATCH_THRESHOLD,
        "mae_threshold": MAE_THRESHOLD,
        "dominance_margin": DOMINANCE_MARGIN,
        "observed_dominance": dominance,
        "sum_exact_match_fraction": candidates["sum"]["exact_match_fraction"],
        "mean_exact_match_fraction": candidates["mean"]["exact_match_fraction"],
        "sum_MAE": candidates["sum"]["mae"],
        "mean_MAE": candidates["mean"]["mae"],
        "sum_median_absolute_error": candidates["sum"]["median_absolute_error"],
        "mean_median_absolute_error": candidates["mean"]["median_absolute_error"],
        "sum_maximum_absolute_error": candidates["sum"]["maximum_absolute_error"],
        "mean_maximum_absolute_error": candidates["mean"]["maximum_absolute_error"],
        "evaluable_documented_windows": winning["evaluable_windows"],
    }


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_input = pd.read_csv(args.cluster_input)
    required_cluster_columns = {"stay_id", "time", "groupHPD", "urineoutput"}
    missing_cluster_columns = required_cluster_columns - set(cluster_input.columns)
    if missing_cluster_columns:
        raise ValueError(f"Missing cluster columns: {sorted(missing_cluster_columns)}")
    if cluster_input.duplicated(["stay_id", "time"]).any():
        raise ValueError("Cluster input contains duplicate patient-time rows")
    unexpected_times = set(cluster_input["time"].dropna().astype(int)) - set(PLANNED_WINDOWS)
    if unexpected_times:
        raise ValueError(f"Unexpected trajectory windows: {sorted(unexpected_times)}")
    included_ids = set(cluster_input["stay_id"].astype(int).unique())
    onsets = pd.read_csv(args.onsets, usecols=["stay_id", "saki_onset"])
    onsets = onsets[onsets["stay_id"].isin(included_ids)].drop_duplicates("stay_id")

    documented = collect_documented_windows(args.events, included_ids)
    documented["charttime_hours"] = documented["charttime"] / 60.0
    documented = documented.merge(onsets, on="stay_id", how="inner", validate="many_to_one")
    raw_window = np.floor((documented["charttime_hours"] - documented["saki_onset"]) / 6.0)
    documented["time"] = np.where(raw_window >= 0, raw_window + 1, raw_window).astype(int)
    documented = documented[documented["time"].between(-2, 28)]

    mask = (
        documented.groupby(["stay_id", "time"], as_index=False)
        .agg(
            uo_documented=("urineoutput", "size"),
            observed_uo_sum=("urineoutput", "sum"),
            observed_uo_mean=("urineoutput", "mean"),
        )
    )
    mask["uo_documented"] = mask["uo_documented"].gt(0)

    audited = cluster_input.merge(mask, on=["stay_id", "time"], how="left", validate="one_to_one")
    audited["uo_documented"] = audited["uo_documented"].eq(True)
    audited["documented_zero"] = audited["uo_documented"] & audited["observed_uo_sum"].eq(0)

    documented_rows = audited["uo_documented"]
    archived_uo = pd.to_numeric(audited.loc[documented_rows, "urineoutput"], errors="coerce")
    aggregation_audit = audit_aggregation_rule(
        archived_uo,
        audited.loc[documented_rows, "observed_uo_sum"],
        audited.loc[documented_rows, "observed_uo_mean"],
    )
    if not aggregation_audit["passed"]:
        failure_status = {
            **aggregation_audit,
            "analysis_inputs_written": False,
            "required_action": (
                "Stop: neither raw-record sum nor mean reconstructs the archived "
                "six-hour urine-output field within the prespecified thresholds."
            ),
        }
        (output_dir / "eicu_uo_scenario_status.json").write_text(
            json.dumps(failure_status, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8",
        )
        print(json.dumps(failure_status, indent=2, ensure_ascii=False))
        return 1

    selected_column = f"observed_uo_{aggregation_audit['selected_rule']}"
    observed_uo = pd.to_numeric(
        audited.loc[documented_rows, selected_column], errors="coerce"
    )
    # Use the selected raw reconstruction explicitly in both sensitivity inputs.
    audited.loc[documented_rows, "urineoutput"] = observed_uo.to_numpy()

    coverage = calculate_coverage(audited)

    documented_only = audited[audited["uo_documented"]].copy()
    eligible = coverage.loc[coverage["documented_windows"] >= 4, "stay_id"]
    documented_only = documented_only[documented_only["stay_id"].isin(eligible)]
    documented_only = documented_only[cluster_input.columns]

    high_coverage_ids = coverage.loc[
        coverage["documented_fraction_planned"] >= 0.50, "stay_id"
    ]
    high_coverage = documented_only[documented_only["stay_id"].isin(high_coverage_ids)].copy()

    mask.to_csv(output_dir / "eicu_uo_documentation_mask.csv", index=False)
    coverage.to_csv(output_dir / "eicu_uo_patient_coverage.csv", index=False)
    documented_only.to_csv(output_dir / "eicu_mixak_documented_windows.csv", index=False)
    high_coverage.to_csv(output_dir / "eicu_mixak_high_coverage.csv", index=False)

    phenotype_names = {1: "DR", 2: "RR", 3: "PW"}
    summary = (
        coverage.assign(phenotype=coverage["groupHPD"].map(phenotype_names))
        .groupby("phenotype")["documented_fraction_planned"]
        .agg(["count", "mean", "median", "min", "max"])
        .reset_index()
    )
    summary.to_csv(output_dir / "eicu_uo_coverage_by_phenotype.csv", index=False)

    prior_available_denominator_n = int(
        coverage["documented_fraction_available"].ge(0.50).sum()
    )
    status = {
        **aggregation_audit,
        "analysis_inputs_written": True,
        "planned_window_definition": list(PLANNED_WINDOWS),
        "planned_window_denominator": len(PLANNED_WINDOWS),
        "authoritative_patients": len(coverage),
        "documented_window_scenario_patients": int(documented_only["stay_id"].nunique()),
        "high_coverage_threshold": ">=15 of 30 planned windows",
        "high_coverage_scenario_patients": int(high_coverage["stay_id"].nunique()),
        "patients_meeting_prior_available_row_fraction": prior_available_denominator_n,
        "patients_excluded_by_corrected_denominator": (
            prior_available_denominator_n - int(high_coverage["stay_id"].nunique())
        ),
        "archived_vs_raw_documented_uo_mismatch_rows": aggregation_audit[
            "discrepancy_windows"
        ],
        "urine_output_aggregation": (
            f"{aggregation_audit['selected_rule']} of documented raw urine-output "
            "records per 6-hour window"
        ),
    }
    (output_dir / "eicu_uo_scenario_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    print(summary.to_string(index=False))
    print(f"Documented-window scenario: {documented_only['stay_id'].nunique()} patients, {len(documented_only)} rows")
    print(f"High-coverage scenario: {high_coverage['stay_id'].nunique()} patients, {len(high_coverage)} rows")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
