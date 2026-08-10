from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
from docx import Document


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "09_revision_documents"
    / "audit_numerical_consistency.py"
)
SPEC = importlib.util.spec_from_file_location("audit_numerical_consistency", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_anchor_tokens_are_required_in_the_same_anchored_paragraph():
    document = Document()
    document.add_paragraph("Primary result: n=123, OR 1.25 (95% CI 1.01–1.55).")
    checks = []
    MODULE.check_anchor_tokens(
        checks,
        check_id="result",
        category="test",
        source="source.csv",
        document_name="document.docx",
        document=document,
        anchor="Primary result",
        tokens=["n=123", "1.25", "1.01–1.55"],
        location="Results",
    )
    assert checks[0].passed is True

    checks = []
    MODULE.check_anchor_tokens(
        checks,
        check_id="result",
        category="test",
        source="source.csv",
        document_name="document.docx",
        document=document,
        anchor="Primary result",
        tokens=["n=123", "0.99"],
        location="Results",
    )
    assert checks[0].passed is False
    assert "0.99" in checks[0].details


def test_repeated_table_headers_select_the_requested_occurrence():
    document = Document()
    for value in ["first", "second"]:
        table = document.add_table(rows=1, cols=2)
        table.rows[0].cells[0].text = "Feature"
        table.rows[0].cells[1].text = "Value"
        cells = table.add_row().cells
        cells[0].text = value
        cells[1].text = "1.000"

    frame = pd.DataFrame([{"Feature": "second", "Value": 1.0}])
    checks = []
    MODULE.check_dataframe_table(
        checks,
        check_id="second_table",
        category="test",
        source="source.csv",
        document_name="document.docx",
        document=document,
        location="second table",
        frame=frame,
        occurrence=1,
    )
    assert checks[0].passed is True
    assert '"matching_tables": 2' in checks[0].details
    assert '"selected_occurrence": 1' in checks[0].details


def test_table_identity_fails_closed_when_occurrence_is_missing():
    document = Document()
    table = document.add_table(rows=1, cols=2)
    table.rows[0].cells[0].text = "Feature"
    table.rows[0].cells[1].text = "Value"
    frame = pd.DataFrame([{"Feature": "row", "Value": 1.0}])
    checks = []
    MODULE.check_dataframe_table(
        checks,
        check_id="missing_table",
        category="test",
        source="source.csv",
        document_name="document.docx",
        document=document,
        location="missing table",
        frame=frame,
        occurrence=1,
    )
    assert checks[0].passed is False


def test_grid_mismatch_reports_exact_cell_coordinates():
    mismatches = MODULE.grid_mismatches(
        [["A", "B"], ["x", "1.000"]],
        [["A", "B"], ["x", "1.001"]],
    )
    assert mismatches == [
        {"row": 1, "column": 1, "expected": "1.000", "observed": "1.001"}
    ]
