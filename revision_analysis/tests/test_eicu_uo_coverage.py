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
