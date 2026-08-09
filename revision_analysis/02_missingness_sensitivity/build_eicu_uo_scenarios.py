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
    observed_uo = pd.to_numeric(
        audited.loc[documented_rows, "observed_uo_mean"], errors="coerce"
    )
    comparison = np.isclose(archived_uo, observed_uo, rtol=1e-9, atol=1e-9, equal_nan=True)
    archived_vs_observed_mismatch = int((~comparison).sum())
    # Use the raw documented-window mean explicitly in both sensitivity inputs.
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
        "status": "PASS",
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
        "archived_vs_raw_documented_uo_mismatch_rows": archived_vs_observed_mismatch,
        "urine_output_aggregation": "mean of documented raw urine-output records per 6-hour window",
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
