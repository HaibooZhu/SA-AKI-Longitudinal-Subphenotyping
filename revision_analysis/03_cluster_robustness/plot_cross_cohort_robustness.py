#!/usr/bin/env python3
"""Plot the completed cross-cohort K=3 processing-sensitivity matrix."""

from __future__ import annotations

from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"


STYLE_DIR = Path(__file__).resolve().parents[1] / "07_tables_figures"
sys.path.insert(0, str(STYLE_DIR))
from publication_figure_style import (  # noqa: E402
    DOUBLE_COLUMN_IN,
    FONT_LEGEND,
    FONT_TEXT,
    NEUTRAL_DARK,
    apply_publication_style,
    export_figure,
)


COHORTS = [("mimic", "MIMIC-IV"), ("eicu", "eICU-CRD"), ("aumc", "AmsterdamUMCdb")]
SCENARIOS = [
    ("zscore_scaling", "Z-score"),
    ("robust_scaling", "Median/IQR"),
    ("exclude_documented_rrt", "Exclude RRT"),
    ("complete_30_window_followup", "Complete\nfollow-up"),
    ("limited_forward_fill_complete_rows", "Limited fill"),
]


def main() -> None:
    repo = Path(__file__).resolve().parents[2]
    report_dir = repo / "02_revision_outputs/reports/W3_cross_cohort_robustness"
    data = pd.read_csv(report_dir / "cross_cohort_robustness_metrics.csv")
    if len(data) != len(COHORTS) * len(SCENARIOS):
        raise ValueError("Expected a complete 3-cohort × 5-scenario robustness grid")

    matrix = np.empty((len(COHORTS), len(SCENARIOS)))
    records: dict[tuple[int, int], pd.Series] = {}
    for row_index, (cohort, _) in enumerate(COHORTS):
        for column_index, (scenario, _) in enumerate(SCENARIOS):
            selected = data.loc[data.cohort.eq(cohort) & data.scenario.eq(scenario)]
            if len(selected) != 1:
                raise ValueError(f"Missing result for {cohort}/{scenario}")
            record = selected.iloc[0]
            matrix[row_index, column_index] = record.ari
            records[(row_index, column_index)] = record

    apply_publication_style()
    fig, ax = plt.subplots(figsize=(DOUBLE_COLUMN_IN, 3.35))
    fig.subplots_adjust(left=0.18, right=0.87, top=0.96, bottom=0.32)
    norm = colors.Normalize(vmin=0.30, vmax=1.00)
    cmap = mpl.colormaps["Blues"]
    ax.set_xticks(range(len(SCENARIOS)), [label for _, label in SCENARIOS])
    ax.set_yticks(range(len(COHORTS)), [label for _, label in COHORTS])
    ax.tick_params(axis="both", length=0, pad=7)
    ax.set_axisbelow(True)
    ax.set_xticks(np.arange(-0.5, len(SCENARIOS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(COHORTS), 1), minor=True)
    ax.grid(which="minor", color="#EFEFEF", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)

    for (row_index, column_index), record in records.items():
        caution = record.screening_threshold_status == "CAUTION"
        bubble_size = 80 + 310 * float(record.retained_fraction)
        facecolor = cmap(norm(record.ari))
        ax.scatter(
            column_index,
            row_index,
            s=bubble_size,
            color=facecolor,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        text_color = "white" if record.ari >= 0.68 else NEUTRAL_DARK
        ax.text(
            column_index,
            row_index,
            f"{record.ari:.2f}",
            ha="center",
            va="center",
            fontsize=FONT_TEXT,
            fontweight="bold",
            color=text_color,
            zorder=4,
        )
        if caution:
            ax.scatter(
                column_index + 0.28,
                row_index - 0.30,
                marker="x",
                s=23,
                linewidth=1.1,
                color="#8A4740",
                zorder=5,
            )

    ax.set_xlim(-0.5, len(SCENARIOS) - 0.5)
    ax.set_ylim(len(COHORTS) - 0.5, -0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

    scalar = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(scalar, ax=ax, fraction=0.035, pad=0.025)
    colorbar.set_label("Adjusted Rand index")
    colorbar.outline.set_linewidth(0.6)

    retained_fractions = (0.25, 0.50, 1.00)
    size_handles = [
        ax.scatter(
            [],
            [],
            s=80 + 310 * fraction,
            color="#9DBAD3",
            edgecolor="white",
            linewidth=0.8,
        )
        for fraction in retained_fractions
    ]
    status_handles = [
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="#9DBAD3", markeredgecolor="white", markersize=6),
        Line2D([], [], marker="x", linestyle="none", color="#8A4740", markersize=5),
    ]
    fig.legend(
        size_handles,
        ["25%", "50%", "100%"],
        title="Patients retained",
        ncol=3,
        loc="lower left",
        bbox_to_anchor=(0.17, 0.01),
        handletextpad=0.3,
        columnspacing=0.9,
        handleheight=2.2,
        borderpad=0.5,
        fontsize=FONT_LEGEND,
        title_fontsize=FONT_LEGEND,
    )
    fig.legend(
        status_handles,
        ["Pass", "Caution"],
        ncol=2,
        loc="lower right",
        bbox_to_anchor=(0.88, 0.035),
        handletextpad=0.35,
        columnspacing=0.9,
        fontsize=FONT_LEGEND,
    )

    export_figure(fig, report_dir / "cross_cohort_robustness_matrix")
    plt.close(fig)


if __name__ == "__main__":
    main()
