from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "02_missingness_sensitivity"
    / "audit_cross_cohort_cluster_inputs.py"
)
SPEC = importlib.util.spec_from_file_location("audit_cross_cohort_cluster_inputs", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_aumc_is_explicitly_three_variable_without_bun():
    spec = MODULE.COHORTS["AUMC"]
    assert spec.renal_features == ("creatinine", "urineoutput", "crea_divide_basecrea")
    assert "bun" not in spec.renal_features


def test_other_cohorts_include_all_four_renal_features():
    expected = {"bun", "creatinine", "urineoutput", "crea_divide_basecrea"}
    assert set(MODULE.COHORTS["MIMIC-IV"].renal_features) == expected
    assert set(MODULE.COHORTS["eICU-CRD"].renal_features) == expected


def test_duplicate_patient_time_rows_are_rejected():
    frame = pd.DataFrame(
        {
            "stay_id": [1, 1],
            "time": [1, 1],
            "groupHPD": [1, 1],
            "bun": [10.0, 10.0],
            "creatinine": [1.0, 1.0],
            "urineoutput": [100.0, 100.0],
            "crea_divide_basecrea": [1.0, 1.0],
        }
    )
    with pytest.raises(ValueError, match="duplicate patient-time"):
        MODULE.validate_schema(frame, MODULE.COHORTS["MIMIC-IV"])


def test_unexpected_cluster_labels_are_rejected():
    frame = pd.DataFrame(
        {
            "stay_id": [1],
            "time": [1],
            "groupHPD": [4],
            "bun": [10.0],
            "creatinine": [1.0],
            "urineoutput": [100.0],
            "crea_divide_basecrea": [1.0],
        }
    )
    with pytest.raises(ValueError, match="unexpected phenotype"):
        MODULE.validate_schema(frame, MODULE.COHORTS["MIMIC-IV"])


def valid_status_frames():
    inventory = pd.DataFrame(
        {"cohort": ["A", "B"], "patient_time_duplicates": [0, 0]}
    )
    baselines = pd.DataFrame(
        {
            "cohort": ["A", "B"],
            "matrix_patients": [10, 20],
            "matrix_patients_linked_to_in_range_baseline": [10, 20],
            "ratio_exact_after_rounding_fraction": [1.0, 0.9995],
            "ratio_absolute_error_max": [0.0, 0.01],
        }
    )
    return inventory, baselines


def test_status_is_computed_from_integrity_thresholds():
    inventory, baselines = valid_status_frames()
    status = MODULE.determine_audit_status(inventory, baselines)
    assert status["passed"] is True
    assert status["failure_reasons"] == []


@pytest.mark.parametrize(
    ("column", "bad_value"),
    [
        ("matrix_patients_linked_to_in_range_baseline", 9),
        ("ratio_exact_after_rounding_fraction", 0.9),
        ("ratio_absolute_error_max", 0.02),
    ],
)
def test_status_fails_when_baseline_audit_misses_threshold(column, bad_value):
    inventory, baselines = valid_status_frames()
    baselines.loc[0, column] = bad_value
    status = MODULE.determine_audit_status(inventory, baselines)
    assert status["passed"] is False
    assert status["failure_reasons"]


def test_status_fails_on_duplicate_patient_time_rows():
    inventory, baselines = valid_status_frames()
    inventory.loc[0, "patient_time_duplicates"] = 1
    status = MODULE.determine_audit_status(inventory, baselines)
    assert status["passed"] is False
