#!/usr/bin/env python3
"""Figure S15 -- documentation-density restrictions select on survival.

Panel a: 28-day mortality among patients removed versus retained by each
documentation restriction, within each trajectory phenotype.
Panel b: retention rate of each phenotype under each restriction.

The figure shows that the restrictions are not ancillary to the outcome: within
Progressive Worsening they preferentially delete the patients who die, so the
restricted populations cannot test whether that phenotype exists.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "01_revision_analysis/07_tables_figures"))

from nature_revision_figure_style import (  # noqa: E402
    NEUTRAL,
    PHENOTYPE,
    WIDTH_DOUBLE_IN,
    apply_style,
    hairline_grid,
    panel_label,
    save_publication,
)

import matplotlib.pyplot as plt  # noqa: E402

REPORT = REPO / "02_revision_outputs/reports/W20_documentation_selection"
FIGURE_DIR = REPO / "02_revision_outputs/figures_redrawn_20260813"

SHORT = {
    ">=50% urine-output coverage": "≥50% urine-output\ncoverage (eICU-CRD)",
    "complete 30-window follow-up": "Complete 30-window follow-up",
}
COHORT_LABEL = {"MIMIC-IV": "MIMIC-IV", "eICU-CRD": "eICU-CRD", "AUMC": "AmsterdamUMCdb"}


def main() -> None:
    apply_style()
    retention = pd.read_csv(REPORT / "restriction_retention.csv")
    tests = pd.read_csv(REPORT / "restriction_selection_tests.csv")
    data = retention.merge(
        tests[["cohort", "restriction", "phenotype", "odds_ratio_death_given_removed", "fisher_p"]],
        on=["cohort", "restriction", "phenotype"],
        how="left",
    )
    data = data[data["phenotype"].ne("ALL")].copy()
    data["scenario"] = data.apply(
        lambda r: f"{COHORT_LABEL[r['cohort']]}\n{SHORT[r['restriction']]}"
        if r["restriction"] == "complete 30-window follow-up"
        else SHORT[r["restriction"]],
        axis=1,
    )
    order = [
        "≥50% urine-output\ncoverage (eICU-CRD)",
        "MIMIC-IV\nComplete 30-window follow-up",
        "eICU-CRD\nComplete 30-window follow-up",
        "AmsterdamUMCdb\nComplete 30-window follow-up",
    ]
    phenotypes = ["RR", "DR", "PW"]

    fig, axes = plt.subplots(
        2, 1, figsize=(WIDTH_DOUBLE_IN, 4.5), sharex=True, height_ratios=[1.35, 1.0]
    )

    # ---- Panel a: mortality of removed vs retained ----
    ax = axes[0]
    group_width, bar_width = 0.78, 0.11
    xticks, xlabels = [], []
    for gi, scenario in enumerate(order):
        block = data[data["scenario"].eq(scenario)].set_index("phenotype")
        for pi, pheno in enumerate(phenotypes):
            if pheno not in block.index:
                continue
            row = block.loc[pheno]
            centre = gi + (pi - 1) * (group_width / 3)
            ax.bar(
                centre - bar_width / 2,
                row["mortality28_removed"],
                bar_width,
                color=PHENOTYPE[pheno],
                edgecolor="none",
                alpha=0.45,
            )
            ax.bar(
                centre + bar_width / 2,
                row["mortality28_retained"],
                bar_width,
                color=PHENOTYPE[pheno],
                edgecolor="none",
            )
            if pheno == "PW" and row["fisher_p"] < 0.05:
                top = max(row["mortality28_removed"], row["mortality28_retained"])
                ax.text(
                    centre,
                    top + 0.045,
                    f"OR {row['odds_ratio_death_given_removed']:.2f}\nP={row['fisher_p']:.2g}",
                    ha="center",
                    va="bottom",
                    fontsize=5.6,
                    color=PHENOTYPE["PW"],
                    linespacing=1.15,
                )
        xticks.append(gi)
        xlabels.append(scenario)
    ax.set_ylim(0, 0.88)
    ax.set_ylabel("28-day mortality")
    hairline_grid(ax)
    panel_label(ax, "a", dx=-0.055)
    ax.set_title(
        "Patients removed by each restriction (pale) versus retained (solid)",
        fontsize=6.8,
        color=NEUTRAL["mid"],
        pad=3,
    )

    # ---- Panel b: retention rate ----
    ax = axes[1]
    for gi, scenario in enumerate(order):
        block = data[data["scenario"].eq(scenario)].set_index("phenotype")
        for pi, pheno in enumerate(phenotypes):
            if pheno not in block.index:
                continue
            centre = gi + (pi - 1) * (group_width / 3)
            ax.bar(
                centre,
                block.loc[pheno, "retention_rate"],
                bar_width * 1.7,
                color=PHENOTYPE[pheno],
                edgecolor="none",
            )
            ax.text(
                centre,
                block.loc[pheno, "retention_rate"] + 0.012,
                f"{int(block.loc[pheno, 'n_retained'])}",
                ha="center",
                va="bottom",
                fontsize=5.6,
                color=NEUTRAL["line"],
            )
    ax.set_ylim(0, 0.42)
    ax.set_ylabel("Fraction retained")
    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels, fontsize=6.2, linespacing=1.2)
    hairline_grid(ax)
    panel_label(ax, "b", dx=-0.055)

    handles = [
        plt.Rectangle((0, 0), 1, 1, color=PHENOTYPE[p], label=n)
        for p, n in [
            ("RR", "Rapid Recovery"),
            ("DR", "Delayed Recovery"),
            ("PW", "Progressive Worsening"),
        ]
    ]
    axes[0].legend(
        handles=handles, loc="upper left", ncol=3, fontsize=6.2, handlelength=1.2,
        columnspacing=1.1, borderaxespad=0.2,
    )
    fig.subplots_adjust(hspace=0.14)
    written = save_publication(fig, FIGURE_DIR, "Figure_S15_documentation_selection")
    for path in written:
        print("wrote", path.relative_to(REPO))


if __name__ == "__main__":
    main()
