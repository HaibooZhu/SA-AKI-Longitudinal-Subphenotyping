from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "05_diuretic_exploratory"
    / "exploratory_diuretic_audit.py"
)
SPEC = importlib.util.spec_from_file_location("exploratory_diuretic_audit", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_landmark_response_definitions_are_explicit():
    full = pd.DataFrame(
        {
            "stay_id": [1, 2, 3, 4],
            "first_use_time": [1, 1, 1, 2],
            "survival_28day": [5, 5, 5, 5],
            "groupHPD": [1, 2, 3, 1],
            "gender": ["M", "F", "M", "F"],
            "age": [60, 61, 62, 63],
            "baseline_Scr": [1, 1, 1, 1],
            "mortality_28d": [0, 0, 1, 0],
        }
    )
    events = pd.DataFrame(
        {
            "stay_id": [1, 2, 3, 4],
            "diuretic_time": [1, 1, 1, 1],
            "diuretic_amout": [20, 20, 20, 20],
            "urineoutput_before_useDiu": [300, 100, 200, 100],
            "urineoutput_after_useDiu": [210, 205, 205, 250],
            "one_label_diu_res": [
                "responsive",
                "responsive",
                "Non-responsive",
                "responsive",
            ],
        }
    )
    result = MODULE.build_landmark("mimic", full, events).set_index("stay_id")
    assert result.index.tolist() == [1, 2, 3]
    assert result.loc[1, "response_archived_first"] == 1
    assert result.loc[1, "response_absolute_200"] == 1
    assert result.loc[1, "response_strict_10pct_200"] == 0
    assert result.loc[2, "response_strict_10pct_200"] == 1
    assert result.loc[3, "response_archived_first"] == 0


def test_matching_spec_documents_phenotype_specific_settings():
    observed = {(row[0], row[1]): row[3] for row in MODULE.MATCHING_SPECS}
    assert observed[("mimic", "DR")] == "caliper=0.05"
    assert observed[("mimic", "RR")] == "M1=1.5; M2=4"
    assert observed[("aumc", "PW")] == "caliper=0.20"
