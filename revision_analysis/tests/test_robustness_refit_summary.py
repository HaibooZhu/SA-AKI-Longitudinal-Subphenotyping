from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "03_cluster_robustness"
    / "summarize_robustness_refits.py"
)
SPEC = importlib.util.spec_from_file_location("summarize_robustness_refits", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_alignment_recovers_permuted_labels():
    original = pd.Series([1, 1, 2, 2, 3, 3])
    refit = pd.Series([3, 3, 1, 1, 2, 2])
    assert MODULE.align_labels(original, refit) == {1: 2, 2: 3, 3: 1}


def test_prespecified_scenario_grid_has_fifteen_refits():
    assert len(MODULE.COHORTS) * len(MODULE.SCENARIOS) == 15
