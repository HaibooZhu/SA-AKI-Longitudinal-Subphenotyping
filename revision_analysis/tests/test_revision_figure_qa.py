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
