from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "03_cluster_robustness"
    / "build_cross_cohort_scenarios.py"
)
SPEC = importlib.util.spec_from_file_location("build_cross_cohort_scenarios", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_zscore_is_patient_independent_global_scaling():
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0]})
    scaled, parameters = MODULE.scale_features(frame, ["x"], "zscore")
    assert np.isclose(scaled.x.mean(), 0)
    assert np.isclose(scaled.x.std(ddof=0), 1)
    assert parameters["x_center"] == 2.5


def test_robust_scaling_uses_median_and_iqr():
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0, 100.0]})
    scaled, parameters = MODULE.scale_features(frame, ["x"], "robust")
    assert np.isclose(scaled.x.median(), 0)
    assert parameters["x_scale"] == 25.5


def test_full_followup_requires_exactly_all_30_planned_windows():
    complete = list(MODULE.PLANNED_WINDOWS)
    frame = pd.DataFrame(
        {
            "stay_id": [1] * 30 + [2] * 29,
            "time": complete + complete[:-1],
        }
    )
    result = MODULE.full_followup(frame)
    assert result.stay_id.unique().tolist() == [1]


def test_rrt_scenario_excludes_only_documented_positive_patients(tmp_path):
    path = tmp_path / "rrt.csv"
    pd.DataFrame({"stay_id": [1, 2], "is_rrt": [1, 0]}).to_csv(path, index=False)
    frame = pd.DataFrame({"stay_id": [1, 2, 3], "time": [1, 1, 1]})
    result, excluded, unknown = MODULE.exclude_documented_rrt(frame, path)
    assert result.stay_id.tolist() == [2, 3]
    assert excluded == 1
    assert unknown == 1
