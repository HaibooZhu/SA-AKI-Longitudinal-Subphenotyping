#!/usr/bin/env python3
"""Fail-closed final traceability audit for the JTIM point-by-point response."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd
from docx import Document


@dataclass(frozen=True)
class Requirement:
    title: str
    result_tokens: tuple[str, ...]
    limitation_tokens: tuple[str, ...]
    evidence_tokens: tuple[str, ...]


REQUIREMENTS = (
    Requirement("E.1 Comprehensive data-integrity audit", ("zero mismatches", "22 mismatches"), ("reporting artifact", "was not used"), ("Tables S3–S5", "data dictionary")),
    Requirement("E.2 Robustness to data-processing assumptions", ("one of three", "415 patients"), ("material sensitivity", "could not be reconstructed"), ("Figure S2", "S15d–S16", "Figures S10–S11")),
    Requirement("E.3 Phenotype definition versus independent validation", ("adjusted ORs", "95% CIs excluded 1"), ("residual confounding",), ("Table S18", "Figure S13")),
    Requirement("E.4 Diuretic-responsiveness analysis", ("0.592", "1.129"), ("exploratory", "treatment effect"), ("Table S19", "Figure S14")),
    Requirement("E.5 Balanced early-classifier evaluation", ("0.836", "0.769", "0.653"), ("no external incremental value", "prospective"), ("Table S17", "Figure S12")),
    Requirement("E.6 Translational claims and positioning", ("cross-cohort reproducibility",), ("future research requirement",), ()),
    Requirement("R1.1 Repositioning and confounding of diuretic responsiveness", ("residual SMDs", "did not replicate"), ("secondary and exploratory", "unmeasured"), ("Table S19",)),
    Requirement("R1.2 Rebalancing classifier performance claims", ("DR-versus-PW", "macro internal/external"), ("weaker",), ("Table S17", "Figure S12")),
    Requirement("R1.3 Critical data-integrity flag in Table S3", ("22 sequential values", "zero source-to-generated mismatches"), ("did not affect",), ("Table S3", "provenance")),
    Requirement("R1.4 Missing urine output and fluid balance", ("1,043 patients", "415 patients", "failed all three"), ("materially limit",), ("Table S16", "Figure S11")),
    Requirement("R1.m1 Baseline creatinine clarification", ("was not imputed", "<0.5 or ≥1.5"), ("selection boundary",), ("Figure S1",)),
    Requirement("R1.m2 Literature-search validity", ("9 August 2026", "2026 SA-AKI"), ("novelty claim was therefore removed",), ("References 31–35",)),
    Requirement("R1.m3 Citation formatting consistency", ("sequential numbered style",), (), ("References",)),
    Requirement("R1.m4 Nonstandard clustering metric", ("2/3 MIMIC-IV", "1/3 eICU-CRD", "2/3 AmsterdamUMCdb"), ("not a validated universal criterion", "does not establish k=3"), ("Table S1", "Table S15e", "Figure S2")),
    Requirement("R1.m5 AmsterdamUMCdb hospital outcome limitation", ("unavailable in AmsterdamUMCdb",), ("restrict cross-cohort outcome interpretation",), ("Table 1",)),
)


def parse_args() -> argparse.Namespace:
    analysis = Path(__file__).resolve().parents[2]
    project = analysis.parents[1]
    package = project / "01_manuscript/14、Journal of Translational Internal Medicine/revision_package_20260805"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--draft", type=Path, default=package / "JTIM_Point-by-point_Response_draft.md")
    parser.add_argument("--response-docx", type=Path, default=package / "JTIM_Point-by-point_Response.docx")
    parser.add_argument("--output-dir", type=Path, default=analysis / "02_revision_outputs/reports/W10_reviewer_traceability")
    parser.add_argument("--package-report", type=Path, default=package / "JTIM_终审_Editor与Reviewer逐条核对_20260810.md")
    return parser.parse_args()


def parse_sections(markdown: str) -> dict[str, str]:
    matches = list(re.finditer(r"^## (.+)$", markdown, flags=re.MULTILINE))
    sections: dict[str, str] = {}
    for index, match in enumerate(matches):
        end = matches[index + 1].start() if index + 1 < len(matches) else len(markdown)
        sections[match.group(1).strip()] = markdown[match.end() : end].strip()
    return sections


def contains_all(text: str, tokens: tuple[str, ...]) -> tuple[bool, list[str]]:
    folded = text.casefold()
    missing = [token for token in tokens if token.casefold() not in folded]
    return not missing, missing


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    sections = parse_sections(args.draft.read_text(encoding="utf-8"))
    word_text = "\n".join(paragraph.text for paragraph in Document(args.response_docx).paragraphs)
    checks: list[dict[str, object]] = []

    for requirement in REQUIREMENTS:
        section = sections.get(requirement.title, "")
        result_ok, missing_results = contains_all(section, requirement.result_tokens)
        limitation_ok, missing_limitations = contains_all(section, requirement.limitation_tokens)
        evidence_ok, missing_evidence = contains_all(section, requirement.evidence_tokens)
        structural = {
            "section_present": bool(section),
            "original_comment_present": "> " in section,
            "response_present": "**Response:**" in section,
            "changes_location_present": "**Changes in the manuscript:**" in section,
            "word_section_present": requirement.title in word_text,
        }
        passed = all(structural.values()) and result_ok and limitation_ok and evidence_ok
        checks.append(
            {
                "comment": requirement.title,
                "status": "PASS" if passed else "FAIL",
                **structural,
                "result_present": result_ok,
                "limitation_present": limitation_ok,
                "supplement_or_manuscript_evidence_present": evidence_ok,
                "missing_result_tokens": "; ".join(missing_results),
                "missing_limitation_tokens": "; ".join(missing_limitations),
                "missing_evidence_tokens": "; ".join(missing_evidence),
            }
        )

    expected_titles = {requirement.title for requirement in REQUIREMENTS}
    extra_sections = sorted(set(sections) - expected_titles)
    missing_sections = sorted(expected_titles - set(sections))
    failures = [check for check in checks if check["status"] == "FAIL"]
    status = {
        "overall_status": "PASS_FINAL_REVIEWER_TRACEABILITY" if not failures and not missing_sections else "FAIL_FINAL_REVIEWER_TRACEABILITY",
        "comments_expected": len(REQUIREMENTS),
        "comments_passed": len(REQUIREMENTS) - len(failures),
        "comments_failed": len(failures),
        "failed_comments": [check["comment"] for check in failures],
        "missing_sections": missing_sections,
        "extra_markdown_sections": extra_sections,
        "submission_readiness": "needs_author_input",
        "author_actions": [
            "approve the final scientific wording, including all transparently reported negative or cautionary findings",
            "confirm author list, affiliations, corresponding-author details, and the response-letter signatory",
            "confirm journal-system metadata and upload-file designations",
            "add final page/line numbers only if the journal specifically requires them after pagination is frozen",
        ],
    }
    pd.DataFrame(checks).to_csv(args.output_dir / "reviewer_traceability_checks.csv", index=False)
    (args.output_dir / "reviewer_traceability_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# JTIM 终审：Editor 与 Reviewer 逐条核对",
        "",
        "版本日期：2026-08-10",
        "",
        f"**自动核对状态：{status['overall_status']}（{status['comments_passed']}/{status['comments_expected']}）**",
        "",
        "本轮按审稿人视角检查每条意见是否同时具备：英文原意见、明确回复、实际结果、必要限制、修改位置，以及对应的正文/补充材料证据。这里的 PASS 表示回复链条完整，不表示所有敏感性分析均得到正向结果。失败、退化或不稳定结果均按原样保留并在回复中解释。",
        "",
        "| 意见 | 原文 | 回复 | 结果 | 限制 | 修改位置/证据 | 终审 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for check in checks:
        yes = lambda value: "是" if value else "否"
        lines.append(
            f"| {check['comment']} | {yes(check['original_comment_present'])} | {yes(check['response_present'])} | {yes(check['result_present'])} | {yes(check['limitation_present'])} | {yes(check['changes_location_present'] and check['supplement_or_manuscript_evidence_present'])} | {check['status']} |"
        )
    lines.extend(
        [
            "",
            "## 综合判断",
            "",
            "- Scientific analyses: complete",
            "- Reviewer-requested robustness analyses: complete",
            "- Repository/evidence integrity: complete",
            "- Final submission QA: complete at the analytical-document level",
            "- Submission readiness: needs_author_input",
            "",
            "当前没有新的统计、方法学或数据完整性阻断项，也不建议为了追求全 PASS 继续增加随机种子、替换模型或追加敏感性实验。eICU 的 K=3 初值敏感性、高尿量覆盖子集 PW 塌缩，以及复杂分类器缺少外部增量价值，均应作为真实研究发现或限制透明保留。",
            "",
            "## 投稿前作者确认",
            "",
        ]
    )
    lines.extend(f"- {action}" for action in status["author_actions"])
    lines.extend(
        [
            "",
            "## QA 证据",
            "",
            "- W8 文档 QA：Table S1 78/78 单元格与 source-of-truth 一致；Figure S2 内嵌资源与发布图 SHA-256 一致。",
            "- W9 数值一致性：54/54 项通过，覆盖正文、标色稿、补充材料、回复信、关键表格和发布图。",
            "- W10 逐条回复追踪：本报告所列 Editor E.1–E.6 与 Reviewer R1.1–R1.m5。",
            "",
        ]
    )
    report = "\n".join(lines)
    (args.output_dir / "W10_FINAL_REVIEWER_TRACEABILITY.md").write_text(report, encoding="utf-8")
    args.package_report.write_text(report, encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))
    return 1 if failures or missing_sections else 0


if __name__ == "__main__":
    raise SystemExit(main())
