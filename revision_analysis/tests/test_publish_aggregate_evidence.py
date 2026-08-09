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
