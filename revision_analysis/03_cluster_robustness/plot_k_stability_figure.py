#!/usr/bin/env python3
"""Plot traceable cross-cohort screening and deep K=2/K=3 evidence for Figure S2."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COHORTS = (("mimic", "MIMIC-IV"), ("eicu", "eICU-CRD"), ("aumc", "AmsterdamUMCdb"))
COLORS = {2: "#3B6FB6", 3: "#C7773E", 4: "#8C8C8C", 5: "#B8B8B8"}


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

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 7.5,
            "axes.titlesize": 9,
            "axes.labelsize": 8,
            "xtick.labelsize": 7,
            "ytick.labelsize": 7,
            "legend.fontsize": 6.5,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 5.2), constrained_layout=True)
    rng = np.random.default_rng(20260809)

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
                s=22,
                color=COLORS[k],
                edgecolor="white",
                linewidth=0.45,
                zorder=3,
            )
            upper.plot(
                [k - 0.16, k + 0.16],
                [np.median(values), np.median(values)],
                color="#202020",
                linewidth=1.1,
                zorder=4,
            )
            selected = int(
                fresh_summary.loc[
                    fresh_summary.cohort.eq(cohort) & fresh_summary.K.eq(k),
                    "selections_across_three_seeds",
                ].iloc[0]
            )
            if selected:
                upper.text(k, max(values) + 0.10, f"selected {selected}/3", ha="center", va="bottom", fontsize=6.2)
        upper.set_title(label, fontweight="bold")
        upper.set_xticks([2, 3, 4, 5])
        upper.set_xlabel("Candidate components (K)")
        upper.set_ylabel("Screening composite score\n(lower is preferred)" if column == 0 else "")
        upper.grid(axis="y", color="#E5E5E5", linewidth=0.6)
        upper.text(-0.17, 1.08, chr(ord("a") + column), transform=upper.transAxes, fontweight="bold", fontsize=9)

        lower = axes[1, column]
        d_cohort = deep.loc[deep.cohort.eq(cohort)].copy()
        for seed, seed_frame in d_cohort.groupby("seed"):
            seed_frame = seed_frame.sort_values("K")
            lower.plot(
                seed_frame.K,
                seed_frame.deep_selection_score,
                color="#9A9A9A",
                linewidth=0.8,
                alpha=0.8,
                zorder=1,
            )
            for row in seed_frame.itertuples():
                lower.scatter(
                    row.K,
                    row.deep_selection_score,
                    s=28,
                    color=COLORS[int(row.K)],
                    edgecolor="white",
                    linewidth=0.5,
                    zorder=3,
                )
        status = str(deep_summary.loc[cohort, "k3_cross_seed_status"])
        status_label = (
            "K=3 reproduced in all 3 starts"
            if status == "PASS_ALL_3_INITIALIZATIONS"
            else "Initialization sensitivity remains"
        )
        lower.text(
            0.5,
            0.96,
            status_label,
            transform=lower.transAxes,
            ha="center",
            va="top",
            fontsize=6.4,
            color="#2E5D34" if status.startswith("PASS") else "#8B3A3A",
        )
        lower.set_xticks([2, 3])
        lower.set_xlim(1.75, 3.25)
        lower.set_xlabel("Candidate components (K)")
        lower.set_ylabel("Deep-fit composite score\n(lower is preferred)" if column == 0 else "")
        lower.grid(axis="y", color="#E5E5E5", linewidth=0.6)
        lower.text(-0.17, 1.08, chr(ord("d") + column), transform=lower.transAxes, fontweight="bold", fontsize=9)

    fig.suptitle(
        "Traceable screening K=2–5 grid and prespecified deep K=2 versus K=3 stability",
        fontsize=10,
        fontweight="bold",
    )
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output_stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(args.output_stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(args.output_stem.with_suffix(".png"), dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(args.output_stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(args.output_stem)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
