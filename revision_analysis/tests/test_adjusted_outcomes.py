from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "04_independent_outcomes"
    / "adjusted_outcomes.py"
)
SPEC = importlib.util.spec_from_file_location("adjusted_outcomes", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_nonrenal_sofa_definition_excludes_renal_component(tmp_path):
    sofa_dir = tmp_path / "04.other_feature_in_three_dataset" / "03.sofa_feature"
    sofa_dir.mkdir(parents=True)
    pd.DataFrame(
        {
            "stay_id": [1, 1, 2],
            "time": [-1, 1, 1],
            "respiration_sofa": [0, 1, 2],
            "coagulation_sofa": [0, 2, 1],
            "liver_sofa": [0, 3, 0],
            "cardiovascular_sofa": [0, 4, 2],
            "cns_sofa": [0, 1, 3],
            "renal_sofa": [0, 4, 4],
        }
    ).to_csv(sofa_dir / "mimic_sofa_clean.csv", index=False)

    result = MODULE.load_onset_nonrenal_sofa(tmp_path, "mimic").set_index("stay_id")
    assert result.loc[1, "onset_nonrenal_sofa"] == 11
    assert result.loc[2, "onset_nonrenal_sofa"] == 8


def test_landmark_populations_do_not_include_day7_deaths_or_missing_status():
    frame = pd.DataFrame(
        {
            "stay_id": [1, 2, 3, 4],
            "mortality_7d": [0, 1, np.nan, 0],
            "observed_through_time28": [True, True, True, False],
        }
    )
    frame["landmark_day7_survivor"] = frame.mortality_7d.eq(0)
    frame["landmark_full_trajectory"] = (
        frame.landmark_day7_survivor & frame.observed_through_time28
    )
    all_survivors = MODULE.select_population(frame, "day7_all_survivors")
    strict = MODULE.select_population(frame, "day7_full_trajectory")
    assert all_survivors.stay_id.tolist() == [1, 4]
    assert strict.stay_id.tolist() == [1]


def test_day8_28_outcome_is_undefined_for_day7_deaths():
    mortality_28d = pd.Series([0.0, 1.0, 1.0, np.nan])
    mortality_7d = pd.Series([0.0, 0.0, 1.0, 0.0])
    result = mortality_28d.where(mortality_7d.eq(0))
    assert result.iloc[0] == 0
    assert result.iloc[1] == 1
    assert np.isnan(result.iloc[2])
    assert np.isnan(result.iloc[3])


def test_mismatch_count_ignores_missing_pairs():
    left = pd.Series([1.0, 2.0, np.nan, 1.0])
    right = pd.Series([1.0, 3.0, 2.0, np.nan])
    assert MODULE.mismatch_count(left, right) == 1


def test_collapse_unique_source_rejects_conflicting_patient_values():
    frame = pd.DataFrame(
        {"stay_id": [1, 1], "first_aki_stage": [1, 2]}
    )
    with np.testing.assert_raises_regex(ValueError, "conflicting first_aki_stage"):
        MODULE.collapse_unique_source(
            frame, ["stay_id", "first_aki_stage"], "stage.csv"
        )
