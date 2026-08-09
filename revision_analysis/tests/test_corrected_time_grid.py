from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "02_missingness_sensitivity"
    / "build_corrected_time_grid.py"
)
SPEC = importlib.util.spec_from_file_location("build_corrected_time_grid", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def fixtures():
    cohort = pd.DataFrame({"stay_id": [1, 2], "groupHPD": [1, 2]})
    observed = pd.DataFrame(
        {"stay_id": [1, 1, 2], "time": [-1, 1, 1], "value": [10.0, None, 99.0]}
    )
    timing = pd.DataFrame(
        {
            "stay_id": [1, 2],
            "intime": [0.0, 0.0],
            "outtime": [8.0, 200.0],
            "saki_onset": [4.0, 20.0],
            "hospital_expire_flag": [0, 1],
        }
    )
    mask = pd.DataFrame(
        {"stay_id": [1, 2], "time": [-1, 1], "uo_documented": [True, True]}
    )
    return cohort, observed, timing, mask


def test_grid_preserves_patients_and_creates_all_planned_windows():
    grid = MODULE.build_planned_grid(*fixtures())
    assert grid.stay_id.nunique() == 2
    assert len(grid) == 60
    assert not grid.duplicated(["stay_id", "time"]).any()
    assert set(grid.groupby("stay_id").size()) == {30}


def test_grid_distinguishes_documented_missing_and_short_followup():
    grid = MODULE.build_planned_grid(*fixtures())
    patient1 = grid[grid.stay_id.eq(1)].set_index("time")
    assert patient1.loc[-1, "window_state"] == "documented_uo"
    assert patient1.loc[1, "window_state"] == "trajectory_row_uo_missing"
    assert patient1.loc[2, "window_state"] == "post_icu_exit_hospital_survivor"
    assert patient1.loc[28, "window_state"] == "post_icu_exit_hospital_survivor"


def test_duplicate_patient_time_is_rejected():
    cohort, observed, timing, mask = fixtures()
    observed = pd.concat([observed, observed.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="duplicate patient-time"):
        MODULE.build_planned_grid(cohort, observed, timing, mask)


def test_forward_fill_never_crosses_patient_boundaries_or_exit_states():
    cohort, observed, timing, mask = fixtures()
    grid = MODULE.build_planned_grid(cohort, observed, timing, mask)
    filled = MODULE.forward_fill_within_patient(
        grid,
        ["value"],
        {"documented_uo", "trajectory_row_uo_missing", "within_icu_no_trajectory_row"},
    )
    patient1 = filled[filled.stay_id.eq(1)].set_index("time")
    patient2 = filled[filled.stay_id.eq(2)].set_index("time")
    assert patient1.loc[1, "value"] == 10.0
    assert pd.isna(patient1.loc[2, "value"])
    assert patient2.loc[1, "value"] == 99.0


def test_all_missing_patient_remains_in_grid():
    cohort, observed, timing, mask = fixtures()
    observed = observed[observed.stay_id.eq(1)]
    mask = mask[mask.stay_id.eq(1)]
    grid = MODULE.build_planned_grid(cohort, observed, timing, mask)
    patient2 = grid[grid.stay_id.eq(2)]
    assert len(patient2) == 30
    assert not patient2["trajectory_row_present"].any()
