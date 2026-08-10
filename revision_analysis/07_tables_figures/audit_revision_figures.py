#!/usr/bin/env python3
"""Fail-closed export and visual-review audit for revised JTIM figures."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree

from PIL import Image


FIGURES = {
    "Figure S2": "W3_deep_k_stability/Figure_S2_cross_cohort_k_stability",
    "Figure S10": "W3_cross_cohort_robustness/cross_cohort_robustness_matrix",
    "Figure S11a": "W3_cluster_robustness/documented_windows_cluster_sensitivity",
    "Figure S11b": "W3_cluster_robustness/high_coverage_cluster_sensitivity",
    "Figure S12": "W4_classifier_validation/Figure_S12_archived_model_calibration_comparator",
    "Figure S13": "W5_independent_outcomes/W5_adjusted_outcomes_forest",
    "Figure S14": "W6_diuretic_exploratory/W6_early_diuretic_response_forest",
}
FORMATS = ("svg", "pdf", "tiff", "png")
MIN_RASTER_EDGE_PX = 1000
MIN_TIFF_DPI = 590
MAX_WIDTH_IN = 7.60
VISUAL_CRITERIA = (
    "conclusion_visible",
    "no_overlap_or_clipping",
    "readable_at_target_width",
    "color_semantics_consistent",
)


@dataclass(frozen=True)
class ArtifactCheck:
    figure_id: str
    format: str
    path: str
    bytes: int
    sha256: str
    width_px: int | None
    height_px: int | None
    dpi_x: float | None
    dpi_y: float | None
    status: str
    detail: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _svg_width_inches(root: ElementTree.Element) -> float | None:
    width = root.attrib.get("width", "")
    match = re.fullmatch(r"([0-9.]+)(pt|in|px)?", width)
    if not match:
        return None
    value = float(match.group(1))
    unit = match.group(2) or "px"
    return {"in": value, "pt": value / 72.0, "px": value / 96.0}[unit]


def inspect_artifact(figure_id: str, path: Path) -> ArtifactCheck:
    if not path.is_file() or path.stat().st_size == 0:
        return ArtifactCheck(
            figure_id, path.suffix.lstrip("."), str(path), 0, "", None, None,
            None, None, "FAIL", "missing or empty",
        )

    fmt = path.suffix.lstrip(".")
    width_px = height_px = None
    dpi_x = dpi_y = None
    failures: list[str] = []
    details: list[str] = []

    if fmt in {"png", "tiff"}:
        with Image.open(path) as image:
            width_px, height_px = image.size
            dpi = image.info.get("dpi", (0, 0))
            dpi_x, dpi_y = float(dpi[0]), float(dpi[1])
        details.append(f"{width_px}x{height_px}px")
        if min(width_px, height_px) < MIN_RASTER_EDGE_PX:
            failures.append(f"short edge below {MIN_RASTER_EDGE_PX}px")
        if fmt == "tiff":
            details.append(f"{dpi_x:.1f}x{dpi_y:.1f}dpi")
            if min(dpi_x, dpi_y) < MIN_TIFF_DPI:
                failures.append(f"TIFF resolution below {MIN_TIFF_DPI}dpi")
    elif fmt == "svg":
        root = ElementTree.parse(path).getroot()
        width_in = _svg_width_inches(root)
        text_nodes = sum(1 for node in root.iter() if node.tag.endswith("text"))
        details.append(f"{text_nodes} editable text nodes")
        if text_nodes == 0:
            failures.append("SVG contains no editable text")
        if width_in is not None:
            details.append(f"{width_in:.2f}in wide")
            if width_in > MAX_WIDTH_IN:
                failures.append(f"width exceeds {MAX_WIDTH_IN:.2f}in")
    elif fmt == "pdf":
        if not path.read_bytes().startswith(b"%PDF"):
            failures.append("invalid PDF signature")

    return ArtifactCheck(
        figure_id=figure_id,
        format=fmt,
        path=str(path),
        bytes=path.stat().st_size,
        sha256=sha256(path),
        width_px=width_px,
        height_px=height_px,
        dpi_x=dpi_x,
        dpi_y=dpi_y,
        status="FAIL" if failures else "PASS",
        detail="; ".join(failures or details or ["valid non-empty export"]),
    )


def validate_visual_review_row(row: dict[str, str], current_png_sha256: str) -> dict[str, object]:
    reviewed_at = row.get("reviewed_at", "").strip()
    try:
        datetime.fromisoformat(reviewed_at.replace("Z", "+00:00"))
        timestamp_valid = True
    except ValueError:
        timestamp_valid = False
    reviewed_sha256 = row.get("reviewed_png_sha256", "").strip().lower()
    hash_format_valid = re.fullmatch(r"[0-9a-f]{64}", reviewed_sha256) is not None
    hash_match = hash_format_valid and reviewed_sha256 == current_png_sha256.lower()
    criteria_pass = all(row.get(field, "") == "PASS" for field in VISUAL_CRITERIA)
    metadata_complete = bool(row.get("reviewed_by", "").strip()) and timestamp_valid
    effective_status = "PASS" if (
        row.get("status") == "PASS"
        and criteria_pass
        and metadata_complete
        and hash_match
    ) else "FAIL"
    return {
        "status": effective_status,
        "reviewed_png_sha256": reviewed_sha256,
        "current_png_sha256": current_png_sha256,
        "hash_match": hash_match,
        "reviewed_by": row.get("reviewed_by", "").strip(),
        "reviewed_at": reviewed_at,
        "criteria_pass": criteria_pass,
        "metadata_complete": metadata_complete,
    }


def load_visual_review(path: Path, report_root: Path) -> dict[str, dict[str, object]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    observed = {row["figure_id"]: row for row in rows}
    if set(observed) != set(FIGURES):
        missing = sorted(set(FIGURES) - set(observed))
        extra = sorted(set(observed) - set(FIGURES))
        raise ValueError(f"visual-review figure mismatch: missing={missing}, extra={extra}")
    return {
        figure_id: validate_visual_review_row(
            observed[figure_id],
            sha256(report_root / f"{stem}.png")
            if (report_root / f"{stem}.png").is_file()
            else "",
        )
        for figure_id, stem in FIGURES.items()
    }


def write_outputs(
    checks: list[ArtifactCheck], visual: dict[str, dict[str, object]], output_dir: Path
) -> str:
    output_dir.mkdir(parents=True, exist_ok=True)
    overall = "PASS" if (
        all(check.status == "PASS" for check in checks)
        and all(review["status"] == "PASS" for review in visual.values())
    ) else "FAIL"

    csv_path = output_dir / "figure_export_manifest.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(asdict(checks[0])), lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(asdict(check) for check in checks)

    status = {
        "overall_status": overall,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "figure_count": len(FIGURES),
        "artifact_count": len(checks),
        "artifact_pass_count": sum(check.status == "PASS" for check in checks),
        "visual_review": visual,
        "criteria": {
            "formats": list(FORMATS),
            "minimum_raster_short_edge_px": MIN_RASTER_EDGE_PX,
            "minimum_tiff_dpi": MIN_TIFF_DPI,
            "maximum_vector_width_in": MAX_WIDTH_IN,
            "svg_text_must_remain_editable": True,
        },
    }
    (output_dir / "figure_qa_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    lines = [
        "# W11 revision figure QA",
        "",
        f"**Overall status: {overall}**",
        "",
        f"- Figures reviewed: {len(FIGURES)}",
        f"- Export artifacts passing: {status['artifact_pass_count']}/{len(checks)}",
        "- Required exports: editable SVG, PDF, 600-dpi TIFF, preview PNG",
        "- Manual visual review: conclusion visibility, clipping/overlap, target-width readability, and color semantics",
        "- Fail-closed binding: every manual PASS is tied to the current PNG SHA-256, reviewer, and review timestamp",
        "",
        "| Figure | Visual review | Export QA |",
        "|---|---:|---:|",
    ]
    for figure_id in FIGURES:
        export_status = "PASS" if all(
            check.status == "PASS" for check in checks if check.figure_id == figure_id
        ) else "FAIL"
        lines.append(f"| {figure_id} | {visual[figure_id]['status']} | {export_status} |")
    lines.extend(
        [
            "",
            "The SHA-256 and technical metadata for every export are recorded in `figure_export_manifest.csv`.",
        ]
    )
    (output_dir / "W11_REVISION_FIGURE_QA.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )
    return overall


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-root", type=Path, required=True)
    parser.add_argument("--visual-review", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    visual = load_visual_review(args.visual_review, args.report_root)
    checks = [
        inspect_artifact(figure_id, args.report_root / f"{stem}.{fmt}")
        for figure_id, stem in FIGURES.items()
        for fmt in FORMATS
    ]
    overall = write_outputs(checks, visual, args.output_dir)
    print(f"W11 revision figure QA: {overall}")
    if overall != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
