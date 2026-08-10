from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1] / "03_cluster_robustness"


def load_module(name: str):
    path = ROOT / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


DEEP_K = load_module("summarize_deep_k_stability")
DEEP_ROBUSTNESS = load_module("summarize_deep_robustness")
UO = load_module("summarize_uo_multiseed_sensitivity")


def test_deep_k_grid_has_eighteen_prespecified_fits():
    assert len(DEEP_K.expected_grid()) == 18


def test_deep_k_alignment_recovers_permuted_labels():
    reference = pd.Series([1, 1, 2, 2, 3, 3])
    candidate = pd.Series([3, 3, 1, 1, 2, 2])
    assert DEEP_K.align_labels(reference, candidate, 3) == {1: 2, 2: 3, 3: 1}


def test_deep_k_score_is_normalized_within_each_cohort_and_seed():
    diagnostics = pd.DataFrame(
        {
            "cohort": ["mimic"] * 4,
            "seed": [1, 1, 2, 2],
            "K": [2, 3, 2, 3],
            "mean_deviance": [100.0, 90.0, 1000.0, 1100.0],
            "high_absolute_lag1_fraction": [0.0, 0.0, 0.0, 0.5],
            "sorted_cluster_prevalence": ["0.5;0.5", "0.3;0.3;0.4"] * 2,
        }
    )
    scored = DEEP_K.score_fits(diagnostics).set_index(["seed", "K"])
    assert scored.loc[(1, 2), "scaled_deviance"] == pytest.approx(1.0)
    assert scored.loc[(1, 3), "scaled_deviance"] == pytest.approx(0.0)
    assert scored.loc[(2, 2), "scaled_deviance"] == pytest.approx(0.0)
    assert scored.loc[(2, 3), "scaled_deviance"] == pytest.approx(1.0)
    assert bool(scored.loc[(1, 3), "selected_within_seed"])
    assert bool(scored.loc[(2, 2), "selected_within_seed"])


def test_all_six_caution_scenarios_have_three_deep_seeds():
    assert len(DEEP_ROBUSTNESS.SCENARIOS) == 6
    assert len(DEEP_ROBUSTNESS.SCENARIOS) * len(DEEP_ROBUSTNESS.SEEDS) == 18


def test_uo_multiseed_grid_has_six_prespecified_fits():
    assert UO.SCENARIOS == ("documented_windows", "high_coverage")
    assert len(UO.SCENARIOS) * len(UO.SEEDS) == 6
