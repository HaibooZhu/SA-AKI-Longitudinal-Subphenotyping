from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "07_tables_figures"
    / "publish_aggregate_evidence.py"
)
SPEC = importlib.util.spec_from_file_location("publish_aggregate_evidence", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_csv_with_patient_identifier_is_rejected(tmp_path):
    path = tmp_path / "unsafe.csv"
    pd.DataFrame({"stay_id": [1], "metric": [0.5]}).to_csv(path, index=False)
    with pytest.raises(ValueError, match="prohibited identifier"):
        MODULE.audit_file(path)


def test_aggregate_csv_without_identifier_passes(tmp_path):
    path = tmp_path / "safe.csv"
    pd.DataFrame({"cohort": ["eICU"], "patients": [1417]}).to_csv(
        path, index=False
    )
    assert MODULE.audit_file(path) == {"rows": 1, "columns": 2}


def test_required_editor_evidence_workstreams_are_whitelisted():
    assert "W1_data_integrity/audit_status.json" in MODULE.EVIDENCE_FILES
    assert (
        "W6_diuretic_exploratory/early_first_dose_interaction_tests.csv"
        in MODULE.EVIDENCE_FILES
    )
    assert "W3_deep_k_stability/deep_k_status.json" in MODULE.EVIDENCE_FILES


def test_workstation_paths_are_sanitized(tmp_path):
    report_root = tmp_path / "revision" / "02_revision_outputs" / "reports"
    report_root.mkdir(parents=True)
    source_path = report_root.parents[1] / "00_frozen_inputs" / "input.csv"
    text = f"Source: `{source_path}`; remote: `/home/researcher/private/input.csv`"
    sanitized = MODULE.sanitize_text(text, report_root)
    assert str(report_root.parents[1]) not in sanitized
    assert "/home/researcher" not in sanitized
    assert "<REVISION_REPOSITORY>" in sanitized
