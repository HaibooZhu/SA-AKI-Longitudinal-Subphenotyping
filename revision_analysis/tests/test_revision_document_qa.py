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


def build_supplement(
    path: Path, include_k6: bool = False, eicu_k3_selections: int = 1
) -> None:
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
            cells[2].text = "3"
            cells[3].text = "100"
            cells[4].text = "0"
            cells[5].text = str(eicu_k3_selections if cohort == "eICU-CRD" and k == 3 else (2 if k == 3 else 0))
    document.save(path)


def test_supplement_accepts_exact_traceable_k_grid(tmp_path):
    path = tmp_path / "supplement.docx"
    build_supplement(path)
    result = MODULE.supplement_checks(path)
    assert result["figure_s2_caption_count"] == 1
    assert result["table_s1_has_only_k2_k5"] is True
    assert result["table_s1_has_three_cohorts"] is True
    assert result["table_s1_selection_annotation_present"] is True
    assert result["table_s1_eicu_k3_not_unanimous"] is True


def test_supplement_rejects_k6(tmp_path):
    path = tmp_path / "supplement.docx"
    build_supplement(path, include_k6=True)
    assert MODULE.supplement_checks(path)["table_s1_has_only_k2_k5"] is False


def test_supplement_rejects_unanimous_eicu_k3_annotation(tmp_path):
    path = tmp_path / "supplement.docx"
    build_supplement(path, eicu_k3_selections=3)
    result = MODULE.supplement_checks(path)
    assert result["table_s1_eicu_k3_not_unanimous"] is False


def test_required_content_checks_cover_all_submission_guardrails():
    text = """
    AmsterdamUMCdb used the other three responses because BUN was unavailable.
    eICU-CRD retained initialization sensitivity.
    The ≥50%-coverage sensitivity failed because the PW component collapsed.
    The classifier showed no external incremental advantage over logistic regression.
    The diuretic analysis remains secondary and exploratory.
    Complete 30-window follow-up was cautionary, with material sensitivity to the
    missing-data rule. The screening grid did not select k=3 unanimously.
    """
    assert all(MODULE.required_content_checks(text).values())


def test_required_content_checks_fail_when_a_guardrail_is_missing():
    checks = MODULE.required_content_checks(
        "AmsterdamUMCdb used the other three responses; eICU retained initialization sensitivity."
    )
    assert checks["aumc_uses_three_clustering_variables"] is True
    assert checks["high_coverage_uo_does_not_preserve_pw"] is False
