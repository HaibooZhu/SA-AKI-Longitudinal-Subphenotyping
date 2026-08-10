"""Shared publication style and export helpers for JTIM revision figures."""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl


SINGLE_COLUMN_IN = 3.50
DOUBLE_COLUMN_IN = 7.20

FONT_PANEL = 8.0
FONT_TEXT = 7.0
FONT_TICK = 6.5
FONT_LEGEND = 6.2

LINE_MAIN = 1.2
LINE_AUX = 0.7
NEUTRAL_DARK = "#3F3F3F"
NEUTRAL_MID = "#777777"
NEUTRAL_LIGHT = "#D9D9D9"
GRID_LIGHT = "#E8E8E8"

# These phenotype colors are fixed across every revised JTIM figure.
PHENOTYPE_COLORS = {
    "DR": "#3B6FB6",
    "RR": "#6A9F58",
    "PW": "#C7773E",
}
K_COLORS = {
    2: "#40566F",
    3: "#738399",
    4: "#9BA6B2",
    5: "#C2C8CE",
}
STATUS_COLORS = {"PASS": "#3F6B48", "CAUTION": "#8A4740"}


def apply_publication_style() -> None:
    """Apply the shared double-column journal figure style."""
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "font.size": FONT_TEXT,
            "axes.titlesize": FONT_TEXT,
            "axes.labelsize": FONT_TEXT,
            "xtick.labelsize": FONT_TICK,
            "ytick.labelsize": FONT_TICK,
            "legend.fontsize": FONT_LEGEND,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.frameon": False,
            "savefig.facecolor": "white",
        }
    )


def add_panel_label(ax, label: str, x: float = -0.13, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        fontsize=FONT_PANEL,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def export_figure(fig, stem: Path, png_dpi: int = 350) -> list[Path]:
    """Export editable vectors plus high-resolution submission and preview rasters."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    outputs = [
        stem.with_suffix(".svg"),
        stem.with_suffix(".pdf"),
        stem.with_suffix(".tiff"),
        stem.with_suffix(".png"),
    ]
    fig.savefig(outputs[0], bbox_inches="tight")
    fig.savefig(outputs[1], bbox_inches="tight")
    fig.savefig(
        outputs[2],
        dpi=600,
        bbox_inches="tight",
        pil_kwargs={"compression": "tiff_lzw"},
    )
    fig.savefig(outputs[3], dpi=png_dpi, bbox_inches="tight")
    return outputs
