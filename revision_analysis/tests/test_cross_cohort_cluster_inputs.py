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
COHORTS = MODULE.cohort_specs(Path("/authorized/project_snapshot"))


def test_aumc_is_explicitly_three_variable_without_bun():
    spec = COHORTS["AUMC"]
    assert spec.renal_features == ("creatinine", "urineoutput", "crea_divide_basecrea")
    assert "bun" not in spec.renal_features


def test_other_cohorts_include_all_four_renal_features():
    expected = {"bun", "creatinine", "urineoutput", "crea_divide_basecrea"}
    assert set(COHORTS["MIMIC-IV"].renal_features) == expected
    assert set(COHORTS["eICU-CRD"].renal_features) == expected


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
        MODULE.validate_schema(frame, COHORTS["MIMIC-IV"])


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
        MODULE.validate_schema(frame, COHORTS["MIMIC-IV"])
