from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

from docx import Document


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "09_revision_documents"
    / "verify_revision_documents.py"
)
SPEC = importlib.util.spec_from_file_location("verify_revision_documents", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def build_supplement(path: Path, include_k6: bool = False) -> None:
    document = Document()
    document.add_paragraph("Figure S2. Traceable candidate-K evidence.")
    table = document.add_table(rows=1, cols=6)
    for index, value in enumerate(
        ["Dataset", "Clusters (K)", "Seeds", "Deviance", "Lag", "Selections"]
    ):
        table.rows[0].cells[index].text = value
    for cohort in ["MIMIC-IV", "eICU-CRD", "AmsterdamUMCdb"]:
        for k in [2, 3, 4, 6 if include_k6 else 5]:
            cells = table.add_row().cells
            cells[0].text = cohort
            cells[1].text = str(k)
    document.save(path)


def test_supplement_accepts_exact_traceable_k_grid(tmp_path):
    path = tmp_path / "supplement.docx"
    build_supplement(path)
    result = MODULE.supplement_checks(path)
    assert result["figure_s2_caption_count"] == 1
    assert result["table_s1_has_only_k2_k5"] is True
    assert result["table_s1_has_three_cohorts"] is True


def test_supplement_rejects_k6(tmp_path):
    path = tmp_path / "supplement.docx"
    build_supplement(path, include_k6=True)
    assert MODULE.supplement_checks(path)["table_s1_has_only_k2_k5"] is False
