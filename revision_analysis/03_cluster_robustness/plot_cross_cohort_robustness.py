#!/usr/bin/env python3
"""Plot the completed cross-cohort K=3 processing-sensitivity matrix."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd


COHORTS = [("mimic", "MIMIC-IV"), ("eicu", "eICU-CRD"), ("aumc", "AmsterdamUMCdb")]
SCENARIOS = [
    ("zscore_scaling", "Z-score\nscaling"),
    ("robust_scaling", "Median/IQR\nscaling"),
    ("exclude_documented_rrt", "Exclude documented\nRRT recipients"),
    ("complete_30_window_followup", "Complete 30-window\nfollow-up"),
    ("limited_forward_fill_complete_rows", "Limited forward fill +\ncomplete renal rows"),
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

    fig, ax = plt.subplots(figsize=(12.4, 4.8))
    fig.subplots_adjust(left=0.14, right=0.90, top=0.84, bottom=0.27)
    norm = colors.Normalize(vmin=0.30, vmax=1.00)
    image = ax.imshow(matrix, cmap="Blues", norm=norm, aspect="auto")
    ax.set_xticks(range(len(SCENARIOS)), [label for _, label in SCENARIOS], fontsize=10)
    ax.set_yticks(range(len(COHORTS)), [label for _, label in COHORTS], fontsize=11)
    ax.tick_params(axis="both", length=0, pad=9)
    ax.set_title(
        "Cross-cohort K=3 processing sensitivity",
        fontsize=15,
        fontweight="bold",
        loc="left",
        pad=14,
    )

    for (row_index, column_index), record in records.items():
        caution = record.screening_threshold_status == "CAUTION"
        annotation = (
            f"Agreement {record.exact_agreement:.2f}\n"
            f"ARI {record.ari:.2f}\n"
            f"Retained {record.retained_fraction:.0%}\n"
            f"{record.screening_threshold_status}"
        )
        text_color = "white" if record.ari >= 0.67 else "#102A43"
        ax.text(
            column_index,
            row_index,
            annotation,
            ha="center",
            va="center",
            fontsize=8.4,
            fontweight="bold" if caution else "normal",
            color=text_color,
            linespacing=1.35,
        )
        if caution:
            ax.add_patch(
                Rectangle(
                    (column_index - 0.49, row_index - 0.49),
                    0.98,
                    0.98,
                    fill=False,
                    edgecolor="#C2413A",
                    linewidth=3.0,
                )
            )

    for boundary in np.arange(-0.5, len(SCENARIOS), 1):
        ax.axvline(boundary, color="white", linewidth=1.5)
    for boundary in np.arange(-0.5, len(COHORTS), 1):
        ax.axhline(boundary, color="white", linewidth=1.5)
    ax.set_xlim(-0.5, len(SCENARIOS) - 0.5)
    ax.set_ylim(len(COHORTS) - 0.5, -0.5)
    for spine in ax.spines.values():
        spine.set_visible(False)

    colorbar = fig.colorbar(image, ax=ax, shrink=0.78, pad=0.02)
    colorbar.set_label("Adjusted Rand index", fontsize=10)
    colorbar.ax.tick_params(labelsize=9)
    fig.text(
        0.14,
        0.04,
        "Red border = prespecified CAUTION (agreement <0.75, ARI <0.50, or minimum cluster prevalence <0.03).",
        fontsize=9,
        color="#5A6470",
    )

    for suffix in ["png", "pdf"]:
        fig.savefig(
            report_dir / f"cross_cohort_robustness_matrix.{suffix}",
            dpi=300,
            bbox_inches="tight",
            facecolor="white",
        )
    plt.close(fig)


if __name__ == "__main__":
    main()
