#!/usr/bin/env python3
"""Build eICU clustering inputs that do not treat undocumented urine output as zero."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


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


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_input = pd.read_csv(args.cluster_input)
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
    audited["uo_documented"] = audited["uo_documented"].fillna(False).astype(bool)
    audited["documented_zero"] = audited["uo_documented"] & audited["observed_uo_sum"].eq(0)

    coverage = (
        audited.groupby(["stay_id", "groupHPD"], as_index=False)
        .agg(
            total_windows=("time", "size"),
            documented_windows=("uo_documented", "sum"),
            documented_zero_windows=("documented_zero", "sum"),
        )
    )
    coverage["documented_fraction"] = coverage["documented_windows"] / coverage["total_windows"]

    documented_only = audited[audited["uo_documented"]].copy()
    eligible = coverage.loc[coverage["documented_windows"] >= 4, "stay_id"]
    documented_only = documented_only[documented_only["stay_id"].isin(eligible)]
    documented_only = documented_only[cluster_input.columns]

    high_coverage_ids = coverage.loc[coverage["documented_fraction"] >= 0.50, "stay_id"]
    high_coverage = documented_only[documented_only["stay_id"].isin(high_coverage_ids)].copy()

    mask.to_csv(output_dir / "eicu_uo_documentation_mask.csv", index=False)
    coverage.to_csv(output_dir / "eicu_uo_patient_coverage.csv", index=False)
    documented_only.to_csv(output_dir / "eicu_mixak_documented_windows.csv", index=False)
    high_coverage.to_csv(output_dir / "eicu_mixak_high_coverage.csv", index=False)

    phenotype_names = {1: "DR", 2: "RR", 3: "PW"}
    summary = (
        coverage.assign(phenotype=coverage["groupHPD"].map(phenotype_names))
        .groupby("phenotype")["documented_fraction"]
        .agg(["count", "mean", "median", "min", "max"])
        .reset_index()
    )
    summary.to_csv(output_dir / "eicu_uo_coverage_by_phenotype.csv", index=False)

    print(summary.to_string(index=False))
    print(f"Documented-window scenario: {documented_only['stay_id'].nunique()} patients, {len(documented_only)} rows")
    print(f"High-coverage scenario: {high_coverage['stay_id'].nunique()} patients, {len(high_coverage)} rows")


if __name__ == "__main__":
    main()
