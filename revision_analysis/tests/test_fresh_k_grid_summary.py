from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "03_cluster_robustness"
    / "summarize_fresh_k_grid.py"
)
SPEC = importlib.util.spec_from_file_location("summarize_fresh_k_grid", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_expected_grid_contains_36_unique_fits():
    assert len(MODULE.expected_grid()) == 36


def test_scoring_is_within_each_cohort_and_seed():
    frame = pd.DataFrame(
        {
            "cohort": ["a"] * 4,
            "seed": [1] * 4,
            "K": [2, 3, 4, 5],
            "mean_deviance": [4.0, 2.0, 3.0, 5.0],
            "high_absolute_lag1_fraction": [0.4, 0.0, 0.2, 0.8],
        }
    )
    result = MODULE.score_grid(frame)
    selected = result.loc[result.selected_within_seed, "K"].tolist()
    assert selected == [3]
