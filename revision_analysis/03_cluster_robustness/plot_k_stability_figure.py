#!/usr/bin/env python3
"""Plot traceable cross-cohort screening and deep K=2/K=3 evidence for Figure S2."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


STYLE_DIR = Path(__file__).resolve().parents[1] / "07_tables_figures"
sys.path.insert(0, str(STYLE_DIR))
from publication_figure_style import (  # noqa: E402
    DOUBLE_COLUMN_IN,
    FONT_LEGEND,
    K_COLORS,
    LINE_AUX,
    LINE_MAIN,
    STATUS_COLORS,
    add_panel_label,
    apply_publication_style,
    export_figure,
)


COHORTS = (("mimic", "MIMIC-IV"), ("eicu", "eICU-CRD"), ("aumc", "AmsterdamUMCdb"))


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--fresh-diagnostics",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_fresh_k_grid/fresh_k_grid_all_seed_diagnostics.csv",
    )
    parser.add_argument(
        "--fresh-summary",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_fresh_k_grid/fresh_k_grid_summary.csv",
    )
    parser.add_argument(
        "--deep-diagnostics",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_deep_k_stability/deep_k_fit_diagnostics.csv",
    )
    parser.add_argument(
        "--deep-summary",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_deep_k_stability/deep_k_cohort_summary.csv",
    )
    parser.add_argument(
        "--output-stem",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_deep_k_stability/Figure_S2_cross_cohort_k_stability",
    )
    return parser.parse_args()


def validate_grid(frame: pd.DataFrame, ks: set[int]) -> None:
    expected = {(cohort, k, seed) for cohort, _ in COHORTS for k in ks for seed in (20260805, 20260806, 20260807)}
    observed = set(zip(frame.cohort, frame.K.astype(int), frame.seed.astype(int)))
    if observed != expected:
        raise ValueError(f"Incomplete K grid; missing={sorted(expected - observed)}")


def main() -> int:
    args = parse_args()
    fresh = pd.read_csv(args.fresh_diagnostics)
    fresh_summary = pd.read_csv(args.fresh_summary)
    deep = pd.read_csv(args.deep_diagnostics)
    deep_summary = pd.read_csv(args.deep_summary).set_index("cohort")
    validate_grid(fresh, {2, 3, 4, 5})
    validate_grid(deep, {2, 3})

    apply_publication_style()
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(DOUBLE_COLUMN_IN, 5.0),
        sharey="row",
        constrained_layout=True,
        gridspec_kw={"height_ratios": [1.0, 0.92]},
    )
    rng = np.random.default_rng(20260809)
    upper_limit = max(1.48, float(fresh.historical_selection_score.max()) + 0.10)
    lower_limit = max(1.48, float(deep.deep_selection_score.max()) + 0.10)

    for column, (cohort, label) in enumerate(COHORTS):
        upper = axes[0, column]
        f_cohort = fresh.loc[fresh.cohort.eq(cohort)]
        for k in (2, 3, 4, 5):
            values = f_cohort.loc[
                f_cohort.K.eq(k), "historical_selection_score"
            ].to_numpy(float)
            jitter = rng.normal(0, 0.035, len(values))
            upper.scatter(
                np.full(len(values), k) + jitter,
                values,
                s=19,
                color=K_COLORS[k],
                edgecolor="white",
                linewidth=0.45,
                zorder=3,
            )
            upper.plot(
                [k - 0.16, k + 0.16],
                [np.median(values), np.median(values)],
                color="#202020",
                linewidth=LINE_MAIN,
                zorder=4,
            )
            selected = int(
                fresh_summary.loc[
                    fresh_summary.cohort.eq(cohort) & fresh_summary.K.eq(k),
                    "selections_across_three_seeds",
                ].iloc[0]
            )
            if selected:
                upper.text(
                    k,
                    max(values) + 0.06,
                    f"{selected}/3",
                    ha="center",
                    va="bottom",
                    fontsize=FONT_LEGEND,
                    color="#555555",
                )
        upper.set_title(label, fontweight="bold", pad=5)
        upper.set_xticks([2, 3, 4, 5])
        upper.set_xlabel("Components, K")
        upper.set_ylabel("Screening composite score\n(lower is preferred)" if column == 0 else "")
        upper.set_ylim(-0.06, upper_limit)
        upper.grid(axis="y", color="#E8E8E8", linewidth=LINE_AUX)
        add_panel_label(upper, chr(ord("a") + column))

        lower = axes[1, column]
        d_cohort = deep.loc[deep.cohort.eq(cohort)].copy()
        for seed, seed_frame in d_cohort.groupby("seed"):
            seed_frame = seed_frame.sort_values("K")
            lower.plot(
                seed_frame.K,
                seed_frame.deep_selection_score,
                color="#9A9A9A",
                linewidth=LINE_AUX,
                alpha=0.8,
                zorder=1,
            )
            for row in seed_frame.itertuples():
                lower.scatter(
                    row.K,
                    row.deep_selection_score,
                    s=24,
                    color=K_COLORS[int(row.K)],
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=3,
                )
        status = str(deep_summary.loc[cohort, "k3_cross_seed_status"])
        stable = status == "PASS_ALL_3_INITIALIZATIONS"
        status_label = "Stable" if stable else "Sensitive"
        lower.text(
            0.03,
            0.95,
            status_label,
            transform=lower.transAxes,
            ha="left",
            va="top",
            fontsize=FONT_LEGEND,
            fontweight="bold",
            color=STATUS_COLORS["PASS" if stable else "CAUTION"],
        )
        lower.set_xticks([2, 3])
        lower.set_xlim(1.75, 3.25)
        lower.set_xlabel("Components, K")
        lower.set_ylabel("Deep-fit composite score\n(lower is preferred)" if column == 0 else "")
        lower.set_ylim(-0.06, lower_limit)
        lower.grid(axis="y", color="#E8E8E8", linewidth=LINE_AUX)
        add_panel_label(lower, chr(ord("d") + column))

    handles = [
        Line2D([], [], marker="o", linestyle="none", color=K_COLORS[k], label=f"K={k}")
        for k in (2, 3, 4, 5)
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.015),
        ncol=4,
        columnspacing=1.0,
        handletextpad=0.3,
    )
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    export_figure(fig, args.output_stem)
    plt.close(fig)
    print(args.output_stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
