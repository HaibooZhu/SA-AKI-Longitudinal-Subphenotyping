#!/usr/bin/env python3
"""Summarize prespecified deep K=2 versus K=3 multi-initialization refits."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


COHORTS = {
    "mimic": "01.MIMICIV_SAKI_trajCluster/df_mixAK_fea4_C3.csv",
    "eicu": "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv",
    "aumc": "02.AUMCdb_SAKI_trajCluster/df_mixAK_fea3_C3_aumc.csv",
}
KS = (2, 3)
SEEDS = (20260805, 20260806, 20260807)
AGREEMENT_THRESHOLD = 0.75
ARI_THRESHOLD = 0.50
MINIMUM_CLUSTER_PREVALENCE = 0.03
MAXIMUM_HIGH_LAG1_FRACTION = 0.0


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--refit-dir",
        type=Path,
        default=repo / "02_revision_outputs/intermediate/W3_deep_k_refits",
    )
    parser.add_argument(
        "--snapshot-root",
        type=Path,
        default=repo / "00_frozen_inputs/data_snapshot/remote_project_snapshot",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W3_deep_k_stability",
    )
    return parser.parse_args()


def expected_grid() -> set[tuple[str, int, int]]:
    return {
        (cohort, k, seed)
        for cohort in COHORTS
        for k in KS
        for seed in SEEDS
    }


def align_labels(reference: pd.Series, candidate: pd.Series, k: int) -> dict[int, int]:
    contingency = pd.crosstab(candidate, reference).reindex(
        index=range(1, k + 1), columns=range(1, k + 1), fill_value=0
    )
    rows, columns = linear_sum_assignment(-contingency.to_numpy())
    return {
        int(contingency.index[row]): int(contingency.columns[column])
        for row, column in zip(rows, columns)
    }


def minmax(values: pd.Series) -> pd.Series:
    lower, upper = values.min(), values.max()
    if np.isclose(lower, upper):
        return pd.Series(0.0, index=values.index)
    return (values - lower) / (upper - lower)


def load_grid(refit_dir: Path) -> tuple[pd.DataFrame, dict[tuple[str, int, int], pd.DataFrame]]:
    diagnostics: list[pd.DataFrame] = []
    assignments: dict[tuple[str, int, int], pd.DataFrame] = {}
    missing: list[str] = []
    for cohort, k, seed in sorted(expected_grid()):
        stem = f"{cohort}__primary_deep__K{k}__seed{seed}"
        diagnostic_path = refit_dir / f"{stem}__diagnostics.csv"
        assignment_path = refit_dir / f"{stem}__assignments.csv"
        if not diagnostic_path.exists() or not assignment_path.exists():
            missing.append(stem)
            continue
        diagnostic = pd.read_csv(diagnostic_path)
        if len(diagnostic) != 1:
            raise ValueError(f"{stem}: expected exactly one diagnostic row")
        diagnostics.append(diagnostic)
        assignments[(cohort, k, seed)] = pd.read_csv(assignment_path)
    if missing:
        raise FileNotFoundError("Missing deep K refits: " + ", ".join(missing))
    frame = pd.concat(diagnostics, ignore_index=True)
    if frame.duplicated(["cohort", "K", "seed"]).any():
        raise ValueError("Duplicate deep K diagnostic rows")
    return frame, assignments


def score_fits(diagnostics: pd.DataFrame) -> pd.DataFrame:
    scored = diagnostics.copy()
    # Match the fresh K-grid definition: compare candidate K values only within
    # the same cohort and initialization. Pooling seeds would change the
    # min-max denominator and make the two panels of Figure S2 incomparable.
    score_groups = scored.groupby(["cohort", "seed"], group_keys=False)
    scored["scaled_deviance"] = score_groups["mean_deviance"].transform(minmax)
    scored["scaled_lag1_failure"] = score_groups[
        "high_absolute_lag1_fraction"
    ].transform(minmax)
    scored["deep_selection_score"] = np.sqrt(
        scored.scaled_deviance**2 + scored.scaled_lag1_failure**2
    )
    scored["minimum_cluster_prevalence"] = scored.sorted_cluster_prevalence.map(
        lambda value: min(float(item) for item in str(value).split(";"))
    )
    scored["diagnostically_admissible"] = (
        scored.high_absolute_lag1_fraction.le(MAXIMUM_HIGH_LAG1_FRACTION)
        & scored.minimum_cluster_prevalence.ge(MINIMUM_CLUSTER_PREVALENCE)
    )
    scored["rank_within_cohort_seed"] = scored.groupby(["cohort", "seed"])[
        "deep_selection_score"
    ].rank(method="min")
    scored["selected_within_seed"] = scored.rank_within_cohort_seed.eq(1)
    return scored


def k3_vs_archived(
    assignments: dict[tuple[str, int, int], pd.DataFrame], snapshot_root: Path
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort, relative in COHORTS.items():
        archived = pd.read_csv(snapshot_root / relative, usecols=["stay_id", "groupHPD"])
        archived = archived.drop_duplicates("stay_id")
        for seed in SEEDS:
            candidate = assignments[(cohort, 3, seed)]
            merged = candidate.merge(archived, on="stay_id", validate="one_to_one")
            mapping = align_labels(merged.groupHPD, merged.group_median, 3)
            aligned = merged.group_median.map(mapping)
            prevalence = aligned.value_counts(normalize=True).reindex([1, 2, 3], fill_value=0)
            rows.append(
                {
                    "cohort": cohort,
                    "seed": seed,
                    "patients": len(merged),
                    "label_mapping": json.dumps(mapping, sort_keys=True),
                    "exact_agreement_vs_archived": float(aligned.eq(merged.groupHPD).mean()),
                    "ari_vs_archived": float(adjusted_rand_score(merged.groupHPD, aligned)),
                    "nmi_vs_archived": float(
                        normalized_mutual_info_score(merged.groupHPD, aligned)
                    ),
                    "minimum_cluster_prevalence": float(prevalence.min()),
                }
            )
    return pd.DataFrame(rows)


def pairwise_k3(
    assignments: dict[tuple[str, int, int], pd.DataFrame]
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cohort in COHORTS:
        for left_seed, right_seed in itertools.combinations(SEEDS, 2):
            left = assignments[(cohort, 3, left_seed)][["stay_id", "group_median"]]
            right = assignments[(cohort, 3, right_seed)][["stay_id", "group_median"]]
            merged = left.merge(
                right,
                on="stay_id",
                suffixes=("_left", "_right"),
                validate="one_to_one",
            )
            mapping = align_labels(merged.group_median_left, merged.group_median_right, 3)
            aligned_right = merged.group_median_right.map(mapping)
            rows.append(
                {
                    "cohort": cohort,
                    "left_seed": left_seed,
                    "right_seed": right_seed,
                    "patients": len(merged),
                    "label_mapping": json.dumps(mapping, sort_keys=True),
                    "exact_agreement": float(
                        merged.group_median_left.eq(aligned_right).mean()
                    ),
                    "ari": float(
                        adjusted_rand_score(merged.group_median_left, aligned_right)
                    ),
                    "nmi": float(
                        normalized_mutual_info_score(
                            merged.group_median_left, aligned_right
                        )
                    ),
                }
            )
    return pd.DataFrame(rows)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    try:
        diagnostics, assignments = load_grid(args.refit_dir)
    except FileNotFoundError as exc:
        status = {"overall_status": "FAIL_INCOMPLETE_DEEP_K_GRID", "reason": str(exc)}
        (args.output_dir / "deep_k_status.json").write_text(
            json.dumps(status, indent=2) + "\n", encoding="utf-8"
        )
        print(json.dumps(status, indent=2))
        return 1

    scored = score_fits(diagnostics)
    archived = k3_vs_archived(assignments, args.snapshot_root)
    pairwise = pairwise_k3(assignments)
    merged_k3 = scored.loc[scored.K.eq(3)].merge(archived, on=["cohort", "seed"])
    merged_k3["fit_pass"] = (
        merged_k3.diagnostically_admissible
        & merged_k3.exact_agreement_vs_archived.ge(AGREEMENT_THRESHOLD)
        & merged_k3.ari_vs_archived.ge(ARI_THRESHOLD)
        & merged_k3.minimum_cluster_prevalence_y.ge(MINIMUM_CLUSTER_PREVALENCE)
    )
    cohort_summary = (
        merged_k3.groupby("cohort", as_index=False)
        .agg(
            k3_fits=("seed", "size"),
            k3_fits_passing=("fit_pass", "sum"),
            median_agreement_vs_archived=("exact_agreement_vs_archived", "median"),
            minimum_ari_vs_archived=("ari_vs_archived", "min"),
            minimum_cluster_prevalence=("minimum_cluster_prevalence_y", "min"),
            maximum_high_lag1_fraction=("high_absolute_lag1_fraction", "max"),
        )
        .merge(
            pairwise.groupby("cohort", as_index=False).agg(
                minimum_pairwise_agreement=("exact_agreement", "min"),
                minimum_pairwise_ari=("ari", "min"),
            ),
            on="cohort",
            validate="one_to_one",
        )
    )
    cohort_summary["k3_cross_seed_status"] = np.where(
        cohort_summary.k3_fits_passing.eq(3)
        & cohort_summary.minimum_pairwise_agreement.ge(AGREEMENT_THRESHOLD)
        & cohort_summary.minimum_pairwise_ari.ge(ARI_THRESHOLD),
        "PASS_ALL_3_INITIALIZATIONS",
        "CAUTION_INITIALIZATION_SENSITIVITY",
    )
    selection_counts = (
        scored.loc[scored.selected_within_seed]
        .groupby(["cohort", "K"])
        .size()
        .unstack(fill_value=0)
        .rename(columns=lambda value: f"K{int(value)}_selections")
        .reset_index()
    )
    cohort_summary = cohort_summary.merge(selection_counts, on="cohort", how="left")
    all_stable = cohort_summary.k3_cross_seed_status.eq(
        "PASS_ALL_3_INITIALIZATIONS"
    ).all()
    status = {
        "overall_status": (
            "DEEP_K3_REPRODUCIBLE_ALL_COHORTS_NONUNIQUE_K"
            if all_stable
            else "DEEP_K3_INITIALIZATION_SENSITIVITY_REMAINS"
        ),
        "expected_fits": len(expected_grid()),
        "observed_fits": len(scored),
        "seeds": list(SEEDS),
        "candidate_k": list(KS),
        "score_normalization_scope": "within cohort and seed across K=2 and K=3",
        "thresholds": {
            "exact_agreement_minimum": AGREEMENT_THRESHOLD,
            "ari_minimum": ARI_THRESHOLD,
            "minimum_cluster_prevalence": MINIMUM_CLUSTER_PREVALENCE,
            "maximum_high_lag1_fraction": MAXIMUM_HIGH_LAG1_FRACTION,
        },
        "cohort_status": cohort_summary[
            ["cohort", "k3_cross_seed_status"]
        ].to_dict(orient="records"),
        "interpretation_guardrail": (
            "Deep K2/K3 comparisons assess initialization stability and reproducibility; "
            "they do not establish a unique biological taxonomy. All prespecified seeds "
            "remain in the denominator."
        ),
    }
    scored.to_csv(args.output_dir / "deep_k_fit_diagnostics.csv", index=False)
    archived.to_csv(args.output_dir / "deep_k_k3_vs_archived.csv", index=False)
    pairwise.to_csv(args.output_dir / "deep_k_k3_pairwise_stability.csv", index=False)
    cohort_summary.to_csv(args.output_dir / "deep_k_cohort_summary.csv", index=False)
    (args.output_dir / "deep_k_status.json").write_text(
        json.dumps(status, indent=2) + "\n", encoding="utf-8"
    )
    report = f"""# Deep K=2 versus K=3 stability experiment

**Status: {status['overall_status']}**

All 18 prespecified deep fits were retained: three cohorts × K=2/3 × three
independent initializations. Each fit used burn 50, keep 2,000, and thin 50.
Numeric labels were aligned before exact-agreement calculations. Conventional
untreated multi-chain R-hat was not calculated because mixture-label switching
precludes direct pooling.

## Cohort-level K=3 reproducibility

{cohort_summary.round(4).to_markdown(index=False)}

## K=3 seed versus archived assignments

{archived.round(4).to_markdown(index=False)}

## Pairwise K=3 seed stability

{pairwise.round(4).to_markdown(index=False)}

Every prespecified initialization, including poorly mixing or degenerate fits, is
reported. The experiment evaluates whether a three-pattern representation can be
reproduced; it does not establish K=3 as the unique true taxonomy. The deep-fit
composite score uses the same normalization scope as the fresh grid: deviance and
lag-1 failure are min-max scaled within each cohort and initialization before K=2
and K=3 are compared.
"""
    (args.output_dir / "W3_DEEP_K2_K3_STABILITY.md").write_text(
        report, encoding="utf-8"
    )
    print(cohort_summary.to_string(index=False))
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
