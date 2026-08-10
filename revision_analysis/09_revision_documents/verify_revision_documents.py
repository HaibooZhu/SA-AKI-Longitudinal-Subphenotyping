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
    "JTIM_Point-by-point_Response.docx",
)
PROHIBITED_EXACT_CLAIMS = (
    "The 3-cluster solution consistently achieved the lowest deviance and acceptable convergence across all datasets.",
    "actionable bedside classifier",
    "supports clinical deployment",
    "precision fluid therapy",
)
CONTENT_FILES = (
    "JTIM_Revised_Manuscript_clean.docx",
    "JTIM_Revised_Manuscript_highlight.docx",
    "JTIM_Point-by-point_Response.docx",
)
REQUIRED_CONTENT_PATTERNS: dict[str, tuple[str, ...]] = {
    "aumc_uses_three_clustering_variables": (
        r"(?:AmsterdamUMCdb|AUMC).{0,450}(?:other three (?:responses|variables)|three (?:renal )?(?:responses|variables))",
    ),
    "eicu_k3_initialization_sensitivity": (
        r"eICU(?:-CRD)?.{0,260}initialization sensitivity",
    ),
    "high_coverage_uo_does_not_preserve_pw": (
        r"(?:≥50%|50%-coverage).{0,500}(?:PW component collapsed|zero (?:median )?PW (?:prevalence|fraction)|PW (?:class|component).{0,60}collapse)",
    ),
    "classifier_has_no_external_incremental_advantage": (
        r"(?:no external incremental (?:value|advantage)|did not improve.{0,140}externally)",
    ),
    "diuretic_analysis_is_exploratory": (
        r"diuretic.{0,450}(?:secondary.{0,80}exploratory|exploratory)",
    ),
    "complete_followup_and_missing_data_cautions": (
        r"complete(?: 30-window)? follow-up.{0,280}(?:caution|sensitivity)",
        r"(?:missing-data.{0,180}(?:sensitivity|rule|caution|refit)|sensitivity.{0,120}missing-data)",
    ),
    "k3_is_not_unanimous_or_uniquely_optimal": (
        r"(?:did not (?:select|provide unanimous support for) k=3|not select k=3 unanimously|k=3.{0,100}not.{0,100}(?:unique|unanimous)|does not establish k=3 as (?:a )?unique)",
    ),
}
NEGATION_GUARDS = re.compile(
    r"\b(?:remove(?:d)?|not|no longer|without|unsupported|cannot|does not)\b",
    flags=re.IGNORECASE,
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


def author_response_text(path: Path) -> str:
    """Return only author-response paragraphs, excluding quoted reviewer wording."""
    document = Document(path)
    prefix = "Response: "
    return "\n".join(
        paragraph.text[len(prefix):]
        for paragraph in document.paragraphs
        if paragraph.text.startswith(prefix)
    )


def prohibited_claim_hits(text: str) -> list[str]:
    """Find legacy claims while ignoring explicit statements that they were removed."""
    hits: list[str] = []
    sentences = re.split(r"(?<=[.!?])\s+|\n+", text)
    for claim in PROHIBITED_EXACT_CLAIMS:
        for sentence in sentences:
            if claim.lower() in sentence.lower() and not NEGATION_GUARDS.search(sentence):
                hits.append(claim)
                break
    return hits


def required_content_checks(text: str) -> dict[str, bool]:
    normalized = re.sub(r"\s+", " ", text)
    return {
        concept: all(re.search(pattern, normalized, flags=re.IGNORECASE) for pattern in patterns)
        for concept, patterns in REQUIRED_CONTENT_PATTERNS.items()
    }


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
    headers = [cell.text.strip() for cell in table_s1.rows[0].cells]
    selection_column = next(
        (index for index, header in enumerate(headers) if "Selections" in header),
        None,
    )
    selection_annotations: dict[str, int] = {}
    if selection_column is not None:
        for row in table_s1.rows[1:]:
            cohort = row.cells[0].text.strip()
            k_value = row.cells[1].text.strip()
            selection = row.cells[selection_column].text.strip()
            if cohort and k_value and selection:
                selection_annotations[f"{cohort}:K{int(float(k_value))}"] = int(float(selection))
    expected_k = [2, 3, 4, 5] * 3
    expected_datasets = {"MIMIC-IV", "eICU-CRD", "AmsterdamUMCdb"}
    expected_selection_keys = {
        f"{cohort}:K{k}" for cohort in expected_datasets for k in (2, 3, 4, 5)
    }
    selection_annotation_present = (
        set(selection_annotations) == expected_selection_keys
        and all(0 <= value <= 3 for value in selection_annotations.values())
    )
    k3_counts = [
        selection_annotations.get(f"{cohort}:K3") for cohort in expected_datasets
    ]
    return {
        "figure_s2_caption_count": len(figure_s2_captions),
        "figure_s2_captions": figure_s2_captions,
        "table_s1_k_values": k_values,
        "table_s1_dataset_values": dataset_values,
        "table_s1_has_only_k2_k5": sorted(k_values) == sorted(expected_k),
        "table_s1_has_three_cohorts": set(dataset_values) == expected_datasets,
        "table_s1_rows": len(table_s1.rows) - 1,
        "table_s1_selection_annotations": selection_annotations,
        "table_s1_selection_annotation_present": selection_annotation_present,
        "table_s1_eicu_k3_not_unanimous": selection_annotations.get("eICU-CRD:K3") != 3,
        "table_s1_k3_not_unanimous_across_cohorts": not all(value == 3 for value in k3_counts),
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
        if not supplement["table_s1_selection_annotation_present"]:
            failures.append("Table S1 must include a complete selections-per-three-seeds annotation")
        if not supplement["table_s1_eicu_k3_not_unanimous"]:
            failures.append("Table S1 must not mark eICU K=3 as unanimously selected")
        if not supplement["table_s1_k3_not_unanimous_across_cohorts"]:
            failures.append("Table S1 must not imply unanimous K=3 selection across cohorts")

    wording_hits: dict[str, list[str]] = {}
    for name in WORDING_FILES:
        path = args.package_dir / name
        if not path.exists():
            continue
        text = author_response_text(path) if name == "JTIM_Point-by-point_Response.docx" else document_text(path)
        hits = prohibited_claim_hits(text)
        wording_hits[name] = hits
        if hits:
            failures.append(f"{name}: prohibited legacy claim remains")
    checks["prohibited_wording_hits"] = wording_hits

    required_content: dict[str, dict[str, bool]] = {}
    for name in CONTENT_FILES:
        path = args.package_dir / name
        if not path.exists():
            continue
        text = author_response_text(path) if name == "JTIM_Point-by-point_Response.docx" else document_text(path)
        concept_checks = required_content_checks(text)
        required_content[name] = concept_checks
        missing_concepts = [concept for concept, passed in concept_checks.items() if not passed]
        if missing_concepts:
            failures.append(f"{name}: missing required scientific content: {', '.join(missing_concepts)}")
    checks["required_scientific_content"] = required_content

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
- Table S1 includes seed-level selection annotations and does not mark eICU K=3 as unanimous: {checks.get('supplement', {}).get('table_s1_selection_annotation_present', False) and checks.get('supplement', {}).get('table_s1_eicu_k3_not_unanimous', False)}
- Prohibited legacy claims detected: {sum(len(value) for value in wording_hits.values())}
- Required scientific content checks passed: {all(all(values.values()) for values in required_content.values()) if required_content else False}

This fail-closed check covers document structure, author-authored wording, and required
scientific interpretation in the clean/highlighted manuscript and point-by-point
response. It complements, but does not replace, rendered-page visual inspection.
"""
    (args.output_dir / "W8_REVISION_DOCUMENT_QA.md").write_text(report, encoding="utf-8")
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
