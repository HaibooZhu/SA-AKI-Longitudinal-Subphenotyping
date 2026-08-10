from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from xml.etree import ElementTree

import matplotlib.pyplot as plt
from PIL import Image


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "07_tables_figures"
    / "publication_figure_style.py"
)
SPEC = importlib.util.spec_from_file_location("publication_figure_style", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_phenotype_colors_are_complete_and_distinct():
    assert set(MODULE.PHENOTYPE_COLORS) == {"DR", "RR", "PW"}
    assert len(set(MODULE.PHENOTYPE_COLORS.values())) == 3


def test_export_contract_preserves_text_and_high_resolution_tiff(tmp_path):
    MODULE.apply_publication_style()
    fig, ax = plt.subplots(figsize=(MODULE.DOUBLE_COLUMN_IN, 2.0))
    ax.plot([0, 1], [0, 1], color=MODULE.PHENOTYPE_COLORS["DR"])
    ax.set_xlabel("Editable label")
    outputs = MODULE.export_figure(fig, tmp_path / "figure")
    plt.close(fig)

    assert {path.suffix for path in outputs} == {".svg", ".pdf", ".tiff", ".png"}
    assert all(path.stat().st_size > 0 for path in outputs)
    svg = ElementTree.parse(tmp_path / "figure.svg").getroot()
    assert any(node.tag.endswith("text") for node in svg.iter())
    with Image.open(tmp_path / "figure.tiff") as image:
        assert min(image.info["dpi"]) >= 599
