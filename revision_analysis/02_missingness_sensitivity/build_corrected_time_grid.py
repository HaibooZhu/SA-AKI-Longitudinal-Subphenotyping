#!/usr/bin/env python3
"""Build the fixed eICU patient-time grid and classify absent windows.

This is a revision-only implementation. It does not alter or replace the archived
preprocessing code. Patient-level grid files remain local analysis inputs; the report
contains aggregated counts only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = REPO_ROOT / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
PLANNED_WINDOWS = (-2, -1, *range(1, 29))


def window_bounds(time: pd.Series) -> tuple[pd.Series, pd.Series]:
    time = pd.to_numeric(time, errors="raise").astype(int)
    start = np.where(time > 0, (time - 1) * 6, time * 6)
    end = np.where(time > 0, time * 6, (time + 1) * 6)
    return pd.Series(start, index=time.index), pd.Series(end, index=time.index)


def validate_patient_time(frame: pd.DataFrame, name: str) -> None:
    required = {"stay_id", "time"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{name} is missing columns: {sorted(missing)}")
    if frame[["stay_id", "time"]].isna().any().any():
        raise ValueError(f"{name} contains missing patient or time identifiers")
    if frame.duplicated(["stay_id", "time"]).any():
        raise ValueError(f"{name} contains duplicate patient-time rows")


def build_planned_grid(
    cohort: pd.DataFrame,
    observed: pd.DataFrame,
    timing: pd.DataFrame,
    documentation_mask: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if cohort["stay_id"].duplicated().any():
        raise ValueError("Cohort contains duplicate stay_id values")
    validate_patient_time(observed, "Observed trajectory")
    if timing["stay_id"].duplicated().any():
        raise ValueError("Timing table contains duplicate stay_id values")
    timing_required = {"stay_id", "intime", "outtime", "saki_onset", "hospital_expire_flag"}
    missing_timing = timing_required - set(timing.columns)
    if missing_timing:
        raise ValueError(f"Timing table is missing columns: {sorted(missing_timing)}")

    patients = cohort[["stay_id", "groupHPD"]].copy()
    if patients["stay_id"].nunique() != len(patients):
        raise ValueError("Authoritative cohort must have one row per patient")
    patients["_join"] = 1
    windows = pd.DataFrame({"time": PLANNED_WINDOWS, "_join": 1})
    grid = patients.merge(windows, on="_join", how="inner", validate="many_to_many").drop(columns="_join")
    expected_rows = len(patients) * len(PLANNED_WINDOWS)
    if len(grid) != expected_rows:
        raise RuntimeError("Planned grid row count is incorrect")

    observed_columns = [column for column in observed.columns if column not in {"groupHPD"}]
    grid = grid.merge(
        observed[observed_columns],
        on=["stay_id", "time"],
        how="left",
        validate="one_to_one",
        indicator="trajectory_merge",
    )
    grid["trajectory_row_present"] = grid["trajectory_merge"].eq("both")
    grid = grid.drop(columns="trajectory_merge")

    if documentation_mask is not None:
        validate_patient_time(documentation_mask, "Documentation mask")
        mask_columns = [
            column
            for column in ["stay_id", "time", "uo_documented"]
            if column in documentation_mask.columns
        ]
        if "uo_documented" not in mask_columns:
            raise ValueError("Documentation mask is missing uo_documented")
        grid = grid.merge(
            documentation_mask[mask_columns],
            on=["stay_id", "time"],
            how="left",
            validate="one_to_one",
        )
        grid["uo_documented"] = grid["uo_documented"].eq(True)
    else:
        grid["uo_documented"] = False

    timing_one = timing[list(timing_required)].copy()
    grid = grid.merge(timing_one, on="stay_id", how="left", validate="many_to_one")
    if grid[["intime", "outtime", "saki_onset"]].isna().any().any():
        raise ValueError("Incomplete ICU entry/exit timing for authoritative patients")
    grid["icu_entry_relative_hours"] = grid["intime"] - grid["saki_onset"]
    grid["icu_exit_relative_hours"] = grid["outtime"] - grid["saki_onset"]
    grid["window_start_hours"], grid["window_end_hours"] = window_bounds(grid["time"])

    start = grid["window_start_hours"]
    end = grid["window_end_hours"]
    entry = grid["icu_entry_relative_hours"]
    exit_time = grid["icu_exit_relative_hours"]
    nonsurvivor = pd.to_numeric(grid["hospital_expire_flag"], errors="coerce").eq(1)

    state = np.full(len(grid), "within_icu_no_trajectory_row", dtype=object)
    state[end.le(entry)] = "pre_icu_entry"
    state[start.lt(entry) & end.gt(entry)] = "icu_entry_during_window"
    after_exit = start.ge(exit_time)
    exit_during = start.lt(exit_time) & end.gt(exit_time)
    state[after_exit & ~nonsurvivor] = "post_icu_exit_hospital_survivor"
    state[after_exit & nonsurvivor] = "post_icu_exit_hospital_nonsurvivor"
    state[exit_during & ~nonsurvivor] = "icu_exit_during_window_hospital_survivor"
    state[exit_during & nonsurvivor] = "icu_exit_during_window_hospital_nonsurvivor"
    state[grid["trajectory_row_present"]] = "trajectory_row_uo_missing"
    state[grid["trajectory_row_present"] & grid["uo_documented"]] = "documented_uo"
    grid["window_state"] = state

    if grid["stay_id"].nunique() != cohort["stay_id"].nunique():
        raise RuntimeError("Patient identifiers were lost while building the grid")
    if grid.duplicated(["stay_id", "time"]).any():
        raise RuntimeError("Grid construction produced duplicate patient-time rows")
    return grid


def forward_fill_within_patient(
    grid: pd.DataFrame, columns: list[str], eligible_states: set[str]
) -> pd.DataFrame:
    """Forward fill only within patient and only across explicitly eligible states."""
    result = grid.sort_values(["stay_id", "time"]).copy()
    original_ids = result["stay_id"].copy()
    for column in columns:
        filled = result.groupby("stay_id", sort=False)[column].ffill()
        eligible = result["window_state"].isin(eligible_states)
        result.loc[eligible, column] = filled.loc[eligible]
    if not result["stay_id"].equals(original_ids):
        raise RuntimeError("stay_id changed during forward fill")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cohort",
        type=Path,
        default=SNAPSHOT / "03.eICU_SAKI_trajCluster/sk_survival.csv",
    )
    parser.add_argument(
        "--observed",
        type=Path,
        default=SNAPSHOT / "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv",
    )
    parser.add_argument(
        "--onsets",
        type=Path,
        default=(
            SNAPSHOT
            / "00.data_eicu/disease_definition/AKI/eicu_saki_event_time.csv"
        ),
    )
    parser.add_argument(
        "--icu-details",
        type=Path,
        default=SNAPSHOT / "00.data_eicu/feature_data/df_eicu_sk_icudetails.csv",
    )
    parser.add_argument(
        "--documentation-mask",
        type=Path,
        default=(
            REPO_ROOT
            / "02_revision_outputs/analysis_inputs/eicu_uo_sensitivity"
            / "eicu_uo_documentation_mask.csv"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "02_revision_outputs/analysis_inputs/eicu_time_grid",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=REPO_ROOT / "02_revision_outputs/reports/W2_missingness_grid",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    cohort = pd.read_csv(args.cohort, usecols=["stay_id", "groupHPD"]).drop_duplicates("stay_id")
    observed = pd.read_csv(args.observed)
    onsets = pd.read_csv(
        args.onsets, usecols=["stay_id", "intime", "outtime", "saki_onset"]
    )
    details = pd.read_csv(args.icu_details, usecols=["stay_id", "hospital_expire_flag"])
    timing = onsets.merge(details, on="stay_id", how="inner", validate="one_to_one")
    timing = timing[timing["stay_id"].isin(cohort["stay_id"])].copy()
    mask = pd.read_csv(args.documentation_mask)

    grid = build_planned_grid(cohort, observed, timing, mask)
    grid.to_csv(args.output_dir / "eicu_planned_30_window_grid.csv", index=False)

    state_counts = (
        grid.groupby(["groupHPD", "window_state"], as_index=False)
        .size()
        .rename(columns={"size": "window_count"})
    )
    state_counts.to_csv(args.report_dir / "eicu_window_state_counts.csv", index=False)
    patient_summary = (
        grid.assign(
            planned_window=1,
            after_exit=grid["window_state"].str.contains("exit"),
        )
        .groupby(["stay_id", "groupHPD"], as_index=False)
        .agg(
            planned_windows=("planned_window", "sum"),
            trajectory_rows=("trajectory_row_present", "sum"),
            documented_uo_windows=("uo_documented", "sum"),
            exit_related_windows=("after_exit", "sum"),
        )
    )
    patient_summary.to_csv(args.output_dir / "eicu_planned_grid_patient_summary.csv", index=False)

    status = {
        "status": "PASS",
        "patients": int(grid["stay_id"].nunique()),
        "planned_windows_per_patient": len(PLANNED_WINDOWS),
        "grid_rows": len(grid),
        "unique_patient_time_rows": int(grid[["stay_id", "time"]].drop_duplicates().shape[0]),
        "patients_with_all_30_trajectory_rows": int(patient_summary["trajectory_rows"].eq(30).sum()),
        "patients_with_trajectory_ending_before_time_28": int(
            observed.groupby("stay_id")["time"].max().lt(28).sum()
        ),
        "state_counts": {
            key: int(value)
            for key, value in grid["window_state"].value_counts().sort_index().items()
        },
        "interpretation": (
            "Absent planned windows are classified by ICU entry/exit timing and hospital "
            "survival status; they are not silently treated as measured zero or as one "
            "homogeneous missing-data mechanism."
        ),
    }
    (args.report_dir / "eicu_time_grid_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report = f"""# W2 eICU 预定时间网格与不完整随访审计

**状态：{status['status']}**

- 权威患者数：{status['patients']:,}
- 每人预定时间窗：{status['planned_windows_per_patient']}（-2、-1、1 至 28）
- 网格总行数：{status['grid_rows']:,}
- 实际拥有全部 30 个轨迹行的患者：{status['patients_with_all_30_trajectory_rows']:,}
- 轨迹在 time 28 前结束的患者：{status['patients_with_trajectory_ending_before_time_28']:,}

缺失窗口已按 ICU 入科前、ICU 内无轨迹行、出 ICU 所在时间窗、出 ICU 后及住院生存状态分类。这里的住院死亡标志并不提供精确死亡时点，因此报告使用“hospital nonsurvivor”而不把所有出 ICU 后窗口误写成“死亡后窗口”。

该网格仅用于修订敏感性分析；历史/公开预处理保持只读。患者级网格位于本地忽略目录，不进入公开 Git 仓库；可公开报告只包含聚合计数。
"""
    (args.report_dir / "W2_EICU_TIME_GRID.md").write_text(report, encoding="utf-8")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
