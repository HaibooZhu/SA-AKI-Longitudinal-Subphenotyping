#!/usr/bin/env python3
"""Fail-closed structural and wording checks for final JTIM revision documents."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from docx import Document


REQUIRED_FILES = (
    "JTIM_Revised_Manuscript_clean.docx",
    "JTIM_Revised_Manuscript_highlight.docx",
    "JTIM_Revised_Supplementary_material.docx",
    "JTIM_Point-by-point_Response.docx",
    "JTIM_逐点回复_中文工作版.docx",
    "JTIM_Revised_Cover_Letter.docx",
    "JTIM_Revised_Highlights.docx",
    "JTIM_Revised_Title_Page.docx",
)
WORDING_FILES = (
    "JTIM_Revised_Manuscript_clean.docx",
    "JTIM_Revised_Manuscript_highlight.docx",
    "JTIM_Revised_Highlights.docx",
    "JTIM_Revised_Cover_Letter.docx",
)
PROHIBITED_EXACT_CLAIMS = (
    "The 3-cluster solution consistently achieved the lowest deviance and acceptable convergence across all datasets.",
    "actionable bedside classifier",
    "supports clinical deployment",
    "precision fluid therapy",
)


def parse_args() -> argparse.Namespace:
    analysis = Path(__file__).resolve().parents[2]
    project = analysis.parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--package-dir",
        type=Path,
        default=(
            project
            / "01_manuscript/14、Journal of Translational Internal Medicine/revision_package_20260805"
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=analysis / "02_revision_outputs/reports/W8_revision_document_qa",
    )
    return parser.parse_args()


def document_text(path: Path) -> str:
    document = Document(path)
    blocks = [paragraph.text for paragraph in document.paragraphs]
    for table in document.tables:
        blocks.extend(cell.text for row in table.rows for cell in row.cells)
    return "\n".join(blocks)


def supplement_checks(path: Path) -> dict[str, object]:
    document = Document(path)
    figure_s2_captions = [
        paragraph.text.strip()
        for paragraph in document.paragraphs
        if re.match(r"^Figure S2(?:\.|\s)", paragraph.text.strip())
    ]
    table_s1 = document.tables[0]
    k_values: list[int] = []
    for row in table_s1.rows[1:]:
        value = row.cells[1].text.strip()
        if value:
            k_values.append(int(float(value)))
    dataset_values = [row.cells[0].text.strip() for row in table_s1.rows[1:]]
    expected_k = [2, 3, 4, 5] * 3
    expected_datasets = {"MIMIC-IV", "eICU-CRD", "AmsterdamUMCdb"}
    return {
        "figure_s2_caption_count": len(figure_s2_captions),
        "figure_s2_captions": figure_s2_captions,
        "table_s1_k_values": k_values,
        "table_s1_dataset_values": dataset_values,
        "table_s1_has_only_k2_k5": sorted(k_values) == sorted(expected_k),
        "table_s1_has_three_cohorts": set(dataset_values) == expected_datasets,
        "table_s1_rows": len(table_s1.rows) - 1,
    }


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    missing = [name for name in REQUIRED_FILES if not (args.package_dir / name).exists()]
    checks: dict[str, object] = {"missing_required_files": missing}
    failures: list[str] = []
    if missing:
        failures.append("Missing required revision documents")

    supplement_path = args.package_dir / "JTIM_Revised_Supplementary_material.docx"
    if supplement_path.exists():
        supplement = supplement_checks(supplement_path)
        checks["supplement"] = supplement
        if supplement["figure_s2_caption_count"] != 1:
            failures.append("Figure S2 must have exactly one caption")
        if not supplement["table_s1_has_only_k2_k5"]:
            failures.append("Table S1 must contain exactly K=2-5 for all three cohorts")
        if not supplement["table_s1_has_three_cohorts"]:
            failures.append("Table S1 cohort set is incomplete")

    wording_hits: dict[str, list[str]] = {}
    for name in WORDING_FILES:
        path = args.package_dir / name
        if not path.exists():
            continue
        text = document_text(path)
        hits = [claim for claim in PROHIBITED_EXACT_CLAIMS if claim.lower() in text.lower()]
        wording_hits[name] = hits
        if hits:
            failures.append(f"{name}: prohibited legacy claim remains")
    checks["prohibited_wording_hits"] = wording_hits

    status = {
        "overall_status": "PASS_FINAL_REVISION_DOCUMENT_QA" if not failures else "FAIL_FINAL_REVISION_DOCUMENT_QA",
        "failures": failures,
        "checks": checks,
    }
    (args.output_dir / "revision_document_qa_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    report = f"""# Final revision document QA

**Status: {status['overall_status']}**

- Required document files present: {not missing}
- Figure S2 caption count: {checks.get('supplement', {}).get('figure_s2_caption_count', 'not checked')}
- Table S1 restricted to traceable K=2–5 across all three cohorts: {checks.get('supplement', {}).get('table_s1_has_only_k2_k5', False)}
- Prohibited legacy claims detected: {sum(len(value) for value in wording_hits.values())}

This check is structural and text-based. It complements, but does not replace,
rendered-page visual inspection of every final Word document.
"""
    (args.output_dir / "W8_REVISION_DOCUMENT_QA.md").write_text(report, encoding="utf-8")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
