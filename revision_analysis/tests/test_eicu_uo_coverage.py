from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "02_missingness_sensitivity"
    / "build_eicu_uo_scenarios.py"
)
SPEC = importlib.util.spec_from_file_location("build_eicu_uo_scenarios", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_planned_window_definition_has_30_windows_and_no_zero():
    assert len(MODULE.PLANNED_WINDOWS) == 30
    assert 0 not in MODULE.PLANNED_WINDOWS
    assert MODULE.PLANNED_WINDOWS[:2] == (-2, -1)
    assert MODULE.PLANNED_WINDOWS[-1] == 28


def test_coverage_uses_fixed_planned_denominator():
    audited = pd.DataFrame(
        {
            "stay_id": [1] * 4 + [2] * 20,
            "groupHPD": [1] * 4 + [2] * 20,
            "time": list(range(1, 5)) + list(range(1, 21)),
            "uo_documented": [True] * 4 + [True] * 15 + [False] * 5,
            "documented_zero": [False] * 24,
        }
    )
    coverage = MODULE.calculate_coverage(audited).set_index("stay_id")
    assert coverage.loc[1, "documented_fraction_available"] == 1.0
    assert coverage.loc[1, "documented_fraction_planned"] == 4 / 30
    assert coverage.loc[2, "documented_fraction_available"] == 15 / 20
    assert coverage.loc[2, "documented_fraction_planned"] == 0.5


def test_aggregation_audit_selects_mean_from_reconstructed_windows():
    archived = pd.Series([10.0, 20.0, 30.0, 40.0])
    result = MODULE.audit_aggregation_rule(
        archived,
        pd.Series([20.0, 40.0, 60.0, 80.0]),
        pd.Series([10.0, 20.0, 30.0, 40.0]),
    )
    assert result["passed"] is True
    assert result["selected_rule"] == "mean"
    assert result["mean_exact_match_fraction"] == 1.0


def test_aggregation_audit_can_select_sum():
    archived = pd.Series([10.0, 20.0, 30.0, 40.0])
    result = MODULE.audit_aggregation_rule(
        archived,
        pd.Series([10.0, 20.0, 30.0, 40.0]),
        pd.Series([5.0, 10.0, 15.0, 20.0]),
    )
    assert result["passed"] is True
    assert result["selected_rule"] == "sum"


def test_aggregation_audit_fails_when_neither_rule_reconstructs_archive():
    archived = pd.Series([10.0, 20.0, 30.0, 40.0])
    result = MODULE.audit_aggregation_rule(
        archived,
        pd.Series([11.0, 21.0, 31.0, 41.0]),
        pd.Series([9.0, 19.0, 29.0, 39.0]),
    )
    assert result["passed"] is False
    assert result["selected_rule"] is None
    assert result["status"] == "FAIL_UO_AGGREGATION_NOT_RECONSTRUCTED"
