#!/usr/bin/env python3
"""Audit the archived diuretic analysis and run a restricted landmark analysis."""

from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import norm


PHENOTYPE = {1: "DR", 2: "RR", 3: "PW"}
COHORT_LABEL = {"mimic": "MIMIC-IV", "aumc": "AUMC"}
RESPONSE_LEVELS = ["No diuretic", "Non-responsive", "responsive"]
MATCHED_VARIABLES = ["creatinine", "urineoutput", "sofa_norenal", "colloid_bolus"]
AUDIT_VARIABLES = MATCHED_VARIABLES + ["age", "male", "weight", "baseline_Scr"]


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=repo / "results/revision/diuretic_exploratory",
    )
    return parser.parse_args()


def load_files(data_dir: Path, cohort: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    folder = data_dir / ("01.mimic" if cohort == "mimic" else "02.aumc")
    full = pd.read_csv(folder / "df_diuretic_responsitive.csv")
    matched = pd.read_csv(folder / "df_diuretic_responsitive_match.csv")
    events = pd.read_csv(folder / "tmp_df_diuretic_responsitive.csv")
    return full, matched, events


def absolute_smd(left: pd.Series, right: pd.Series) -> float:
    left = pd.to_numeric(left, errors="coerce").dropna()
    right = pd.to_numeric(right, errors="coerce").dropna()
    if len(left) < 2 or len(right) < 2:
        return np.nan
    pooled = np.sqrt((left.var(ddof=1) + right.var(ddof=1)) / 2)
    if pooled == 0:
        return 0.0 if left.mean() == right.mean() else np.inf
    return float(abs(left.mean() - right.mean()) / pooled)


def balance_table(cohort: str, full: pd.DataFrame, matched: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for label, frame in [("before_matching", full), ("after_matching", matched)]:
        frame = frame.copy()
        frame["male"] = frame.gender.map({"M": 1.0, "F": 0.0, 1: 1.0, 0: 0.0})
        for group in [1, 2, 3]:
            subset = frame.loc[pd.to_numeric(frame.groupHPD, errors="coerce") == group]
            for variable in AUDIT_VARIABLES:
                pair_values = []
                for a, b in combinations(RESPONSE_LEVELS, 2):
                    value = absolute_smd(
                        subset.loc[subset.label_diu_res == a, variable],
                        subset.loc[subset.label_diu_res == b, variable],
                    )
                    pair_values.append((a, b, value))
                finite = [item for item in pair_values if np.isfinite(item[2])]
                if finite:
                    worst = max(finite, key=lambda x: x[2])
                    max_smd, worst_pair = worst[2], f"{worst[0]} vs {worst[1]}"
                else:
                    max_smd, worst_pair = np.nan, "not estimable"
                rows.append(
                    {
                        "cohort": cohort,
                        "cohort_label": COHORT_LABEL[cohort],
                        "stage": label,
                        "phenotype": PHENOTYPE[group],
                        "variable": variable,
                        "used_in_archived_matching": variable in MATCHED_VARIABLES,
                        "max_absolute_pairwise_smd": max_smd,
                        "worst_pair": worst_pair,
                    }
                )
    return pd.DataFrame(rows)


def build_landmark(cohort: str, full: pd.DataFrame, events: pd.DataFrame) -> pd.DataFrame:
    first_event = (
        events.sort_values(["stay_id", "diuretic_time"])
        .drop_duplicates("stay_id", keep="first")
        .copy()
    )
    frame = full.merge(first_event, on="stay_id", how="inner", validate="one_to_one")
    frame = frame.loc[
        (pd.to_numeric(frame.first_use_time, errors="coerce") == 1)
        & (pd.to_numeric(frame.survival_28day, errors="coerce") > 1)
        & (pd.to_numeric(frame.urineoutput_before_useDiu, errors="coerce") > 0)
        & (pd.to_numeric(frame.urineoutput_after_useDiu, errors="coerce") > 0)
        & (pd.to_numeric(frame.diuretic_amout, errors="coerce") > 0)
    ].copy()
    frame["response_archived_first"] = (frame.one_label_diu_res == "responsive").astype(int)
    frame["response_strict_10pct_200"] = (
        (frame.urineoutput_after_useDiu > 200)
        & (frame.urineoutput_after_useDiu > 1.10 * frame.urineoutput_before_useDiu)
    ).astype(int)
    frame["response_absolute_200"] = (frame.urineoutput_after_useDiu > 200).astype(int)
    frame["groupHPD"] = pd.to_numeric(frame.groupHPD, errors="coerce").astype(int)
    frame["male"] = frame.gender.map({"M": 1.0, "F": 0.0, 1: 1.0, 0: 0.0})
    frame["age10"] = pd.to_numeric(frame.age, errors="coerce") / 10
    frame["log_baseline_scr"] = np.log(pd.to_numeric(frame.baseline_Scr, errors="coerce"))
    frame["log_pre_uo"] = np.log1p(pd.to_numeric(frame.urineoutput_before_useDiu, errors="coerce"))
    frame["log_dose"] = np.log1p(pd.to_numeric(frame.diuretic_amout, errors="coerce"))
    frame["cohort"] = cohort
    return frame


def fit_pooled_sensitivity(cohort: str, frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for response in [
        "response_archived_first",
        "response_strict_10pct_200",
        "response_absolute_200",
    ]:
        formula = (
            f"mortality_28d ~ {response} + C(groupHPD, Treatment(reference=2)) + "
            "age10 + male + log_baseline_scr + sofa_norenal + log_pre_uo + log_dose"
        )
        required = [
            "mortality_28d", response, "groupHPD", "age10", "male", "log_baseline_scr",
            "sofa_norenal", "log_pre_uo", "log_dose"
        ]
        complete = frame.dropna(subset=required).copy()
        fit = smf.glm(formula, data=complete, family=sm.families.Binomial()).fit(cov_type="HC3")
        low, high = fit.conf_int().loc[response]
        rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABEL[cohort],
                "response_definition": response,
                "adjusted_or_response_vs_nonresponse": np.exp(fit.params[response]),
                "ci_low": np.exp(low),
                "ci_high": np.exp(high),
                "p_value": fit.pvalues[response],
                "n_complete": len(complete),
                "deaths": int(complete.mortality_28d.sum()),
                "responsive_n": int(complete[response].sum()),
            }
        )
    return pd.DataFrame(rows)


def fit_interaction(cohort: str, frame: pd.DataFrame) -> tuple[pd.DataFrame, float]:
    response = "response_archived_first"
    formula = (
        f"mortality_28d ~ {response} * C(groupHPD, Treatment(reference=2)) + "
        "age10 + male + log_baseline_scr + sofa_norenal + log_pre_uo + log_dose"
    )
    required = [
        "mortality_28d", response, "groupHPD", "age10", "male", "log_baseline_scr",
        "sofa_norenal", "log_pre_uo", "log_dose"
    ]
    complete = frame.dropna(subset=required).copy()
    fit = smf.glm(formula, data=complete, family=sm.families.Binomial()).fit(cov_type="HC3")
    names = list(fit.params.index)
    base_term = response
    interaction_terms = {
        group: next(
            name for name in names if name.startswith(f"{response}:") and f"[T.{group}]" in name
        )
        for group in [1, 3]
    }
    rows = []
    for group in [1, 2, 3]:
        contrast = np.zeros(len(names))
        contrast[names.index(base_term)] = 1
        if group in interaction_terms:
            contrast[names.index(interaction_terms[group])] = 1
        beta = float(contrast @ fit.params.to_numpy())
        variance = float(contrast @ fit.cov_params().to_numpy() @ contrast)
        se = np.sqrt(max(variance, 0))
        z = beta / se if se > 0 else np.nan
        rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABEL[cohort],
                "phenotype": PHENOTYPE[group],
                "group": group,
                "adjusted_or_response_vs_nonresponse": np.exp(beta),
                "ci_low": np.exp(beta - 1.96 * se),
                "ci_high": np.exp(beta + 1.96 * se),
                "p_value": 2 * norm.sf(abs(z)),
                "n_complete": len(complete),
                "deaths": int(complete.mortality_28d.sum()),
            }
        )
    constraint = ", ".join(f"{term} = 0" for term in interaction_terms.values())
    interaction_p = float(fit.wald_test(constraint, scalar=True).pvalue)
    return pd.DataFrame(rows), interaction_p


def descriptive_landmark(cohort: str, frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby(["groupHPD", "response_archived_first"])
        .mortality_28d.agg(n="size", deaths="sum", mortality_risk="mean")
        .reset_index()
        .assign(
            cohort=cohort,
            cohort_label=COHORT_LABEL[cohort],
            phenotype=lambda x: x.groupHPD.map(PHENOTYPE),
            response=lambda x: x.response_archived_first.map({0: "nonresponsive", 1: "responsive"}),
        )
    )


def plot_interactions(interactions: pd.DataFrame, out_dir: Path) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
        }
    )
    colors = {"DR": "#3B6FB6", "RR": "#6A9F58", "PW": "#C7773E"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), constrained_layout=True)
    for panel, cohort in enumerate(["mimic", "aumc"]):
        ax = axes[panel]
        subset = interactions.loc[interactions.cohort == cohort]
        for y, phenotype in enumerate(["DR", "RR", "PW"]):
            row = subset.loc[subset.phenotype == phenotype].iloc[0]
            ax.errorbar(
                row.adjusted_or_response_vs_nonresponse,
                y,
                xerr=[
                    [row.adjusted_or_response_vs_nonresponse - row.ci_low],
                    [row.ci_high - row.adjusted_or_response_vs_nonresponse],
                ],
                fmt="o",
                color=colors[phenotype],
                capsize=2,
                markersize=4,
                linewidth=1,
            )
        ax.axvline(1, color="#777777", linestyle="--", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_yticks(range(3), ["DR", "RR", "PW"])
        ax.invert_yaxis()
        ax.set_xlabel("Adjusted OR for 28-day mortality\nresponse vs nonresponse (95% CI)")
        ax.set_title(COHORT_LABEL[cohort])
        ax.grid(axis="x", color="#DDDDDD", linewidth=0.5)
        ax.text(-0.12, 1.06, chr(ord("a") + panel), transform=ax.transAxes, fontweight="bold", fontsize=8)
    stem = out_dir / "W6_early_diuretic_response_forest"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(
    out_dir: Path,
    balance: pd.DataFrame,
    descriptive: pd.DataFrame,
    pooled: pd.DataFrame,
    interactions: pd.DataFrame,
    interaction_tests: pd.DataFrame,
) -> None:
    def md(frame: pd.DataFrame) -> str:
        return frame.round(3).to_markdown(index=False)

    balance_summary = (
        balance.groupby(["cohort_label", "stage", "used_in_archived_matching"])
        .max_absolute_pairwise_smd.agg(max_smd="max", median_smd="median", variables="size")
        .reset_index()
    )
    text = f"""# W6 Diuretic analysis audit and exploratory landmark sensitivity analysis

## Bottom line

The archived diuretic analysis cannot identify a treatment effect. It classifies
patients using urine output after treatment, combines responses across an unrestricted
number of post-onset administrations using an ad hoc majority/tie rule, and compares
those post-treatment groups with patients who did not receive a diuretic. The archived
three-group propensity matching used only day-1 creatinine, urine output, nonrenal SOFA,
and colloid input and did not retain matched-triplet identifiers. The diuretic findings
must therefore be demoted to exploratory, associative supplementary results.

## Archived matching balance audit

Maximum pairwise absolute standardized mean differences (SMDs) were calculated among
the three archived groups within each phenotype, before and after archived matching.
An SMD above 0.10 indicates residual imbalance that should not be dismissed by a
nonsignificant p value.

{md(balance_summary)}

Full variable-level results are in `archived_psm_balance_smd.csv`. Variables not used
in the archived matching include age, sex, weight, and baseline creatinine. Other
important unavailable or unaligned confounders include treatment indication, clinician
assessment of congestion, pre-onset fluid exposure, contemporaneous vasopressor dose,
nephrotoxin exposure, dose route, treatment limitation, and the exact timing of all
covariates relative to the first dose.

## Restricted 24-hour landmark cohort

To reduce, but not eliminate, immortal-time and exposure-timing bias, the sensitivity
analysis included only patients whose first post-SA-AKI diuretic administration occurred
in day 1, who were alive beyond 24 hours, and who had positive recorded urine volumes in
both two-hour windows around the first administration. The comparison is response versus
nonresponse among treated patients; untreated patients are not used as a causal control.

{md(descriptive[["cohort_label", "phenotype", "response", "n", "deaths", "mortality_risk"]])}

## Pooled adjusted response association under alternate response rules

Models adjust for phenotype, age, sex, baseline creatinine, day-1 nonrenal SOFA,
pre-dose two-hour urine output, and first-dose amount. The original first-dose rule is
shown alongside stricter alternatives to expose threshold dependence.

{md(pooled.drop(columns="cohort"))}

## Phenotype-specific exploratory associations

{md(interactions.drop(columns=["cohort", "group"]))}

Joint interaction tests:

{md(interaction_tests)}

## Interpretation and reporting decision

1. Remove diuretic responsiveness from the title, abstract conclusion, and main claim.
2. Do not use the no-diuretic group to estimate benefit or harm.
3. Present the restricted first-dose landmark analysis only in the supplement, labeled
   post-treatment prognostic association.
4. State that a standardized prospective furosemide stress protocol, pre-treatment
   eligibility criteria, treatment indication, and time-varying confounding control
   would be required to evaluate treatment effect or effect modification.
5. Do not recommend phenotype-guided diuretic treatment from these data.
"""
    (out_dir / "W6_DIURETIC_EXPLORATORY.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    balance_parts, landmark_parts, descriptive_parts, pooled_parts, interaction_parts = [], [], [], [], []
    interaction_tests = []
    for cohort in ["mimic", "aumc"]:
        full, matched, events = load_files(args.data_dir, cohort)
        balance_parts.append(balance_table(cohort, full, matched))
        landmark = build_landmark(cohort, full, events)
        landmark_parts.append(landmark)
        descriptive_parts.append(descriptive_landmark(cohort, landmark))
        pooled_parts.append(fit_pooled_sensitivity(cohort, landmark))
        interaction, p_value = fit_interaction(cohort, landmark)
        interaction_parts.append(interaction)
        interaction_tests.append(
            {
                "cohort_label": COHORT_LABEL[cohort],
                "test": "response-by-phenotype interaction (2 df)",
                "p_value": p_value,
            }
        )

    balance = pd.concat(balance_parts, ignore_index=True)
    landmark_all = pd.concat(landmark_parts, ignore_index=True)
    descriptive = pd.concat(descriptive_parts, ignore_index=True)
    pooled = pd.concat(pooled_parts, ignore_index=True)
    interactions = pd.concat(interaction_parts, ignore_index=True)
    interaction_tests_df = pd.DataFrame(interaction_tests)

    balance.to_csv(args.out_dir / "archived_psm_balance_smd.csv", index=False)
    landmark_all.to_csv(args.out_dir / "early_first_dose_landmark_cohort.csv", index=False)
    descriptive.to_csv(args.out_dir / "early_first_dose_descriptive.csv", index=False)
    pooled.to_csv(args.out_dir / "early_first_dose_pooled_models.csv", index=False)
    interactions.to_csv(args.out_dir / "early_first_dose_interaction_models.csv", index=False)
    interaction_tests_df.to_csv(args.out_dir / "early_first_dose_interaction_tests.csv", index=False)
    plot_interactions(interactions, args.out_dir)
    write_report(out_dir=args.out_dir, balance=balance, descriptive=descriptive, pooled=pooled,
                 interactions=interactions, interaction_tests=interaction_tests_df)
    print(pooled.to_string(index=False))


if __name__ == "__main__":
    main()
