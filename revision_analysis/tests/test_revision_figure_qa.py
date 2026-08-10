from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from PIL import Image


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "07_tables_figures"
    / "audit_revision_figures.py"
)
SPEC = importlib.util.spec_from_file_location("audit_revision_figures", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_missing_artifact_fails_closed(tmp_path):
    result = MODULE.inspect_artifact("Figure S2", tmp_path / "missing.svg")
    assert result.status == "FAIL"


def test_raster_resolution_contract(tmp_path):
    path = tmp_path / "small.tiff"
    Image.new("RGB", (500, 500), "white").save(path, dpi=(300, 300))
    result = MODULE.inspect_artifact("Figure S2", path)
    assert result.status == "FAIL"
    assert "below" in result.detail


def test_svg_editable_text_contract(tmp_path):
    path = tmp_path / "editable.svg"
    mpl.rcParams["svg.fonttype"] = "none"
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.text(0.5, 0.5, "text")
    fig.savefig(path)
    plt.close(fig)
    result = MODULE.inspect_artifact("Figure S2", path)
    assert result.status == "PASS"


def test_visual_review_pass_is_bound_to_current_png_hash(tmp_path):
    path = tmp_path / "figure.png"
    Image.new("RGB", (40, 30), "white").save(path)
    digest = MODULE.sha256(path)
    row = {
        "reviewed_png_sha256": digest,
        "reviewed_by": "Codex visual QA",
        "reviewed_at": "2026-08-10T12:00:00-04:00",
        "conclusion_visible": "PASS",
        "no_overlap_or_clipping": "PASS",
        "readable_at_target_width": "PASS",
        "color_semantics_consistent": "PASS",
        "status": "PASS",
    }
    assert MODULE.validate_visual_review_row(row, digest)["status"] == "PASS"


def test_visual_review_fails_when_figure_changes_after_review(tmp_path):
    path = tmp_path / "figure.png"
    Image.new("RGB", (40, 30), "white").save(path)
    reviewed_digest = MODULE.sha256(path)
    Image.new("RGB", (40, 30), "black").save(path)
    row = {
        "reviewed_png_sha256": reviewed_digest,
        "reviewed_by": "Codex visual QA",
        "reviewed_at": "2026-08-10T12:00:00-04:00",
        "conclusion_visible": "PASS",
        "no_overlap_or_clipping": "PASS",
        "readable_at_target_width": "PASS",
        "color_semantics_consistent": "PASS",
        "status": "PASS",
    }
    result = MODULE.validate_visual_review_row(row, MODULE.sha256(path))
    assert result["hash_match"] is False
    assert result["status"] == "FAIL"
