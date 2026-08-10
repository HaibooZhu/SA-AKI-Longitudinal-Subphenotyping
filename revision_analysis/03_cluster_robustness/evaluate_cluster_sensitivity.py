#!/usr/bin/env python3
"""Align and evaluate a mixAK sensitivity clustering against archived labels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"


STYLE_DIR = Path(__file__).resolve().parents[1] / "07_tables_figures"
sys.path.insert(0, str(STYLE_DIR))
from publication_figure_style import (  # noqa: E402
    DOUBLE_COLUMN_IN,
    FONT_LEGEND,
    LINE_AUX,
    LINE_MAIN,
    PHENOTYPE_COLORS,
    add_panel_label,
    apply_publication_style,
    export_figure,
)


PHENOTYPE = {1: "DR", 2: "RR", 3: "PW"}
FEATURE_LABEL = {
    "bun": "BUN (mg/dL)",
    "creatinine": "Creatinine (mg/dL)",
    "urineoutput": "Urine output (mL/6 h)",
    "crea_divide_basecrea": "Creatinine / baseline creatinine",
}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--scenario", required=True)
    parser.add_argument("--scenario-input", type=Path, required=True)
    parser.add_argument(
        "--original-data",
        type=Path,
        default=(
            repo
            / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
            / "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv"
        ),
    )
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_cluster_robustness",
    )
    return parser.parse_args()


def align_labels(original: np.ndarray, sensitivity: np.ndarray) -> dict[int, int]:
    contingency = pd.crosstab(
        pd.Series(sensitivity, name="sensitivity"),
        pd.Series(original, name="original"),
    ).reindex(index=[1, 2, 3], columns=[1, 2, 3], fill_value=0)
    rows, cols = linear_sum_assignment(-contingency.to_numpy())
    return {int(contingency.index[r]): int(contingency.columns[c]) for r, c in zip(rows, cols)}


def plot_results(
    scenario: str,
    merged: pd.DataFrame,
    trajectory: pd.DataFrame,
    result_dir: Path,
) -> None:
    apply_publication_style()
    colors = {group: PHENOTYPE_COLORS[name] for group, name in PHENOTYPE.items()}
    fig, axes = plt.subplots(
        3,
        2,
        # Keep the complete six-panel figure and its caption on one landscape
        # supplement page at the 183-mm publication width.
        figsize=(DOUBLE_COLUMN_IN, 6.45),
        constrained_layout=True,
        gridspec_kw={"height_ratios": [0.92, 1.0, 1.0]},
    )
    ax = axes[0, 0]
    table = pd.crosstab(merged.original_group, merged.aligned_group).reindex(
        index=[1, 2, 3], columns=[1, 2, 3], fill_value=0
    )
    normalized = table.div(table.sum(axis=1), axis=0)
    pale_blues = LinearSegmentedColormap.from_list(
        "jtim_pale_blues", ["#F7F9FC", "#C5DBEC", "#6DA7CF"]
    )
    image = ax.imshow(normalized, vmin=0, vmax=1, cmap=pale_blues)
    for row in range(3):
        for column in range(3):
            color = "white" if normalized.iloc[row, column] > 0.70 else "#303030"
            ax.text(column, row, f"{normalized.iloc[row, column]:.2f}", ha="center", va="center", color=color)
    ax.set_xticks(range(3), [PHENOTYPE[x] for x in [1, 2, 3]])
    ax.set_yticks(range(3), [PHENOTYPE[x] for x in [1, 2, 3]])
    ax.set_xlabel("Sensitivity cluster")
    ax.set_ylabel("Archived phenotype")
    ax.set_title("Row-normalized agreement")
    ax.spines[:].set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, shrink=0.68, fraction=0.045, pad=0.04)
    colorbar.outline.set_linewidth(0.6)

    ax = axes[0, 1]
    for group in [1, 2, 3]:
        values = merged.loc[merged.aligned_group == group, "max_median_probability"]
        ax.hist(values, bins=np.linspace(0, 1, 21), histtype="step", linewidth=LINE_MAIN,
                color=colors[group], label=f"{PHENOTYPE[group]} (n={len(values)})")
    ax.axvline(0.5, color="#777777", linestyle="--", linewidth=LINE_AUX)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("Maximum posterior membership probability")
    ax.set_ylabel("Patients")
    ax.set_title("Assignment confidence")
    ax.legend(fontsize=FONT_LEGEND, loc="upper left")

    features = ["bun", "creatinine", "urineoutput", "crea_divide_basecrea"]
    for position, feature in enumerate(features):
        ax = axes[1 + position // 2, position % 2]
        for group in [1, 2, 3]:
            sensitivity = (
                trajectory.loc[trajectory.aligned_group == group]
                .groupby("time")[feature].mean()
            )
            archived = (
                trajectory.loc[trajectory.original_group == group]
                .groupby("time")[feature].mean()
            )
            ax.plot(sensitivity.index, sensitivity.values, color=colors[group], linewidth=LINE_MAIN,
                    marker="o", markersize=2.4, label=f"{PHENOTYPE[group]} sensitivity")
            ax.plot(archived.index, archived.values, color=colors[group], linewidth=LINE_AUX,
                    linestyle="--", alpha=0.75)
        ax.set_xlabel("6-h window relative to SA-AKI onset")
        ax.set_ylabel(FEATURE_LABEL[feature])
        ax.set_title(FEATURE_LABEL[feature])
    for label, axis in zip(["a", "b", "c", "d", "e", "f"], axes.ravel()):
        add_panel_label(axis, label)
    style_handles = [
        Line2D([], [], color="#666666", lw=LINE_MAIN, marker="o", markersize=3, label="Sensitivity"),
        Line2D([], [], color="#888888", lw=LINE_AUX, ls="--", label="Archived"),
    ]
    fig.legend(
        handles=style_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.012),
        ncol=2,
        columnspacing=1.1,
        handlelength=2.0,
    )
    stem = result_dir / f"{scenario}_cluster_sensitivity"
    export_figure(fig, stem)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    args.result_dir.mkdir(parents=True, exist_ok=True)
    assignments = pd.read_csv(args.result_dir / f"{args.scenario}_assignments.csv")
    diagnostics = pd.read_csv(args.result_dir / f"{args.scenario}_diagnostics.csv")
    if len(diagnostics) != 1:
        raise ValueError("Expected exactly one diagnostic row for the scenario")
    convergence_status = diagnostics.loc[0, "convergence_status"]
    probability_min = float(assignments["max_median_probability"].min())
    probability_max = float(assignments["max_median_probability"].max())
    if probability_min < 0 or probability_max > 1:
        raise ValueError("Posterior assignment probabilities are outside 0-1")
    original_long = pd.read_csv(args.original_data)
    original = original_long[["stay_id", "groupHPD"]].drop_duplicates("stay_id")
    original = original.rename(columns={"groupHPD": "original_group"})
    merged = assignments.merge(original, on="stay_id", how="inner", validate="one_to_one")
    mapping = align_labels(merged.original_group.to_numpy(), merged.group_median.to_numpy())
    merged["aligned_group"] = merged.group_median.map(mapping)
    merged["aligned_hpd_group"] = merged.group_hpd.map(mapping).fillna(4).astype(int)
    merged["exact_match"] = merged.aligned_group == merged.original_group
    merged["confident_match"] = np.where(
        merged.aligned_hpd_group == 4,
        np.nan,
        merged.aligned_hpd_group == merged.original_group,
    )

    all_original_n = original.stay_id.nunique()
    confident = merged.loc[merged.aligned_hpd_group != 4]
    metrics = pd.DataFrame(
        [
            {
                "scenario": args.scenario,
                "patients_in_original_eicu": all_original_n,
                "patients_in_sensitivity": len(merged),
                "retained_fraction": len(merged) / all_original_n,
                "label_mapping": json.dumps(mapping, sort_keys=True),
                "agreement_all": merged.exact_match.mean(),
                "ari_all": adjusted_rand_score(merged.original_group, merged.aligned_group),
                "nmi_all": normalized_mutual_info_score(merged.original_group, merged.aligned_group),
                "uncertain_patients": int((merged.aligned_hpd_group == 4).sum()),
                "uncertain_fraction": float((merged.aligned_hpd_group == 4).mean()),
                "confident_patients": len(confident),
                "agreement_confident": (
                    float((confident.aligned_hpd_group == confident.original_group).mean())
                    if len(confident) else np.nan
                ),
                "ari_confident": (
                    adjusted_rand_score(confident.original_group, confident.aligned_hpd_group)
                    if len(confident) else np.nan
                ),
                "convergence_status": convergence_status,
                "posterior_probability_min": probability_min,
                "posterior_probability_max": probability_max,
            }
        ]
    )
    confusion = pd.crosstab(
        merged.original_group.map(PHENOTYPE),
        merged.aligned_group.map(PHENOTYPE),
        margins=True,
    )
    metrics.to_csv(args.result_dir / f"{args.scenario}_comparison_metrics.csv", index=False)
    confusion.to_csv(args.result_dir / f"{args.scenario}_aligned_confusion.csv")
    merged.to_csv(args.result_dir / f"{args.scenario}_aligned_assignments.csv", index=False)

    scenario_long = pd.read_csv(args.scenario_input)
    trajectory = scenario_long.merge(
        merged[["stay_id", "original_group", "aligned_group"]],
        on="stay_id",
        how="inner",
        validate="many_to_one",
    )
    plot_results(args.scenario, merged, trajectory, args.result_dir)

    diagnostic_interpretation = (
        "The lag-1 diagnostic met the prespecified single-chain threshold. This does "
        "not remove the single-chain or selected-subset limitations."
        if str(convergence_status).upper() == "PASS"
        else "A convergence warning prevents this scenario from being presented as "
        "confirmatory evidence even when agreement metrics appear favorable."
    )

    report = f"""# Cluster sensitivity: {args.scenario}

## Bottom line

This analysis re-fits a three-component multivariate longitudinal mixture using the
archived mixAK settings after changing the urine-output documentation rule. Numeric
mixture labels were aligned to the archived phenotypes by maximum overlap before any
agreement metric was calculated.

**Run diagnostic status: {convergence_status}.** {diagnostic_interpretation}

{metrics.round(3).to_markdown(index=False)}

## Aligned patient counts

{confusion.to_markdown()}

## Interpretation guardrails

- Retention and phenotype-specific urine-output documentation must be considered
  alongside label agreement; a high agreement among a selected subset does not prove
  robustness in all original patients.
- Patients whose 95% HPD lower bound for every assigned component did not exceed 0.5
  are reported as uncertain rather than forced into a phenotype.
- Dashed trajectories in the accompanying figure show archived labels within the same
  retained patients; solid trajectories show the sensitivity assignments.
- This is a single-chain sensitivity refit matching the executable archived settings.
  It does not validate the manuscript's prior three-chain/10,000-iteration description,
  which must be corrected.
- Posterior probabilities are used on their native 0-1 scale; no division by two or
  other display-only rescaling is applied.
"""
    (args.result_dir / f"{args.scenario}_CLUSTER_SENSITIVITY.md").write_text(report, encoding="utf-8")
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
