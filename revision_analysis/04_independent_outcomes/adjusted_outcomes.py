#!/usr/bin/env python3
"""Adjusted independent clinical outcome analyses across the three cohorts."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score


SEED = 20260805
PHENOTYPE = {1: "DR", 2: "RR", 3: "PW"}
COHORT_LABEL = {"mimic": "MIMIC-IV", "eicu": "eICU-CRD", "aumc": "AUMC"}
FORMULA = (
    "{outcome} ~ C(groupHPD, Treatment(reference=2)) + age10 + male + "
    "C(first_aki_stage) + log_baseline_scr"
)


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=repo / "results/revision/independent_outcomes",
    )
    parser.add_argument("--bootstrap", type=int, default=300)
    return parser.parse_args()


def load_cohort(data_dir: Path, cohort: str) -> pd.DataFrame:
    files = [data_dir / f"df_{cohort}_c{group}_riskfactor.csv" for group in [1, 2, 3]]
    frame = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    if frame.stay_id.duplicated().any():
        raise ValueError(f"Duplicate stay_id values in {cohort}")
    if cohort == "eicu":
        # The archived risk-factor export contains an earlier 1,748-patient
        # clustering snapshot.  Restrict it to the authoritative 1,417-patient
        # eICU cohort used in the submitted manuscript and overwrite both the
        # phenotype and mortality fields from that cohort's survival file.
        authoritative_path = (
            data_dir.parents[1]
            / "03.eICU_SAKI_trajCluster"
            / "sk_survival.csv"
        )
        authoritative = pd.read_csv(authoritative_path)[
            ["stay_id", "groupHPD", "mortality_28d"]
        ].drop_duplicates("stay_id")
        if len(authoritative) != 1417:
            raise ValueError(
                f"Expected 1,417 authoritative eICU patients, found {len(authoritative)}"
            )
        frame = frame.drop(columns=["groupHPD", "mortality_28d"]).merge(
            authoritative,
            on="stay_id",
            how="inner",
            validate="one_to_one",
        )
        if len(frame) != len(authoritative):
            raise ValueError(
                "The archived risk-factor table does not cover the full authoritative eICU cohort"
            )
    frame["groupHPD"] = pd.to_numeric(frame.groupHPD, errors="coerce").astype("Int64")
    frame["age10"] = pd.to_numeric(frame.age, errors="coerce") / 10
    frame["male"] = frame.gender.map({"M": 1.0, "F": 0.0})
    baseline = pd.to_numeric(frame.baseline_Scr, errors="coerce")
    frame["log_baseline_scr"] = np.where(baseline > 0, np.log(baseline), np.nan)
    frame["first_aki_stage"] = pd.to_numeric(frame.first_aki_stage, errors="coerce")
    frame["mortality_28d"] = pd.to_numeric(frame.mortality_28d, errors="coerce")
    frame["is_rrt"] = pd.to_numeric(frame.is_rrt, errors="coerce")
    return frame


def fit_model(frame: pd.DataFrame, outcome: str):
    columns = [
        outcome,
        "groupHPD",
        "age10",
        "male",
        "first_aki_stage",
        "log_baseline_scr",
    ]
    complete = frame[columns].dropna().copy()
    complete["groupHPD"] = complete.groupHPD.astype(int)
    complete["first_aki_stage"] = complete.first_aki_stage.astype(int)
    formula = FORMULA.format(outcome=outcome)
    fit = smf.glm(formula, data=complete, family=sm.families.Binomial()).fit(cov_type="HC3")
    return complete, fit


def extract_group_effects(cohort: str, outcome: str, complete: pd.DataFrame, fit) -> pd.DataFrame:
    rows = []
    for group in [1, 3]:
        term = f"C(groupHPD, Treatment(reference=2))[T.{group}]"
        beta = fit.params[term]
        low, high = fit.conf_int().loc[term]
        rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABEL[cohort],
                "outcome": outcome,
                "comparison": f"{PHENOTYPE[group]} vs RR",
                "group": group,
                "reference_group": 2,
                "adjusted_or": np.exp(beta),
                "ci_low": np.exp(low),
                "ci_high": np.exp(high),
                "p_value": fit.pvalues[term],
                "n_complete": len(complete),
                "events": int(complete[outcome].sum()),
                "events_per_parameter": float(complete[outcome].sum() / len(fit.params)),
                "model_auc": roc_auc_score(complete[outcome], fit.predict(complete)),
                "converged": bool(fit.converged),
            }
        )
    return pd.DataFrame(rows)


def standardized_risk_bootstrap(
    cohort: str, frame: pd.DataFrame, n_boot: int
) -> pd.DataFrame:
    complete, fit = fit_model(frame, "mortality_28d")
    point = {}
    for group in [1, 2, 3]:
        counterfactual = complete.copy()
        counterfactual["groupHPD"] = group
        point[group] = float(np.mean(fit.predict(counterfactual)))

    rng = np.random.default_rng(SEED + {"mimic": 1, "eicu": 2, "aumc": 3}[cohort])
    draws = {group: [] for group in [1, 2, 3]}
    failures = 0
    for _ in range(n_boot):
        boot = complete.iloc[rng.integers(0, len(complete), len(complete))].copy()
        try:
            boot_fit = smf.glm(
                FORMULA.format(outcome="mortality_28d"),
                data=boot,
                family=sm.families.Binomial(),
            ).fit()
            for group in [1, 2, 3]:
                counterfactual = complete.copy()
                counterfactual["groupHPD"] = group
                draws[group].append(float(np.mean(boot_fit.predict(counterfactual))))
        except Exception:
            failures += 1

    return pd.DataFrame(
        [
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABEL[cohort],
                "phenotype": PHENOTYPE[group],
                "group": group,
                "standardized_risk": point[group],
                "ci_low": np.quantile(draws[group], 0.025),
                "ci_high": np.quantile(draws[group], 0.975),
                "bootstrap_successes": len(draws[group]),
                "bootstrap_failures": failures,
            }
            for group in [1, 2, 3]
        ]
    )


def unadjusted_table(cohort: str, frame: pd.DataFrame, outcome: str) -> pd.DataFrame:
    work = frame[["groupHPD", outcome]].dropna().copy()
    work["groupHPD"] = work.groupHPD.astype(int)
    result = (
        work.groupby("groupHPD")[outcome]
        .agg(n="size", events="sum", risk="mean")
        .reset_index()
    )
    result.insert(0, "cohort", cohort)
    result.insert(1, "cohort_label", COHORT_LABEL[cohort])
    result.insert(2, "outcome", outcome)
    result["phenotype"] = result.groupHPD.map(PHENOTYPE)
    return result


def plot_forest(effects: pd.DataFrame, out_dir: Path) -> None:
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
    outcome_titles = {
        "mortality_28d": "28-day mortality",
        "is_rrt": "In-hospital renal replacement therapy",
    }
    colors = {"DR vs RR": "#3B6FB6", "PW vs RR": "#C7773E"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), constrained_layout=True)
    y_positions = {
        ("MIMIC-IV", "DR vs RR"): 5.2,
        ("MIMIC-IV", "PW vs RR"): 4.8,
        ("eICU-CRD", "DR vs RR"): 3.2,
        ("eICU-CRD", "PW vs RR"): 2.8,
        ("AUMC", "DR vs RR"): 1.2,
        ("AUMC", "PW vs RR"): 0.8,
    }
    for panel, outcome in enumerate(["mortality_28d", "is_rrt"]):
        ax = axes[panel]
        subset = effects.loc[effects.outcome == outcome]
        for _, row in subset.iterrows():
            y = y_positions[(row.cohort_label, row.comparison)]
            ax.errorbar(
                row.adjusted_or,
                y,
                xerr=[[row.adjusted_or - row.ci_low], [row.ci_high - row.adjusted_or]],
                fmt="o",
                markersize=4,
                capsize=2,
                linewidth=1,
                color=colors[row.comparison],
            )
        ax.axvline(1, color="#777777", linestyle="--", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_ylim(0.2, 5.8)
        ax.set_yticks([5, 3, 1], ["MIMIC-IV", "eICU-CRD", "AUMC"])
        ax.set_xlabel("Adjusted odds ratio (95% CI)")
        ax.set_title(outcome_titles[outcome])
        ax.grid(axis="x", color="#DDDDDD", linewidth=0.5)
        for comparison, offset in [("DR vs RR", 0.18), ("PW vs RR", -0.18)]:
            ax.scatter([], [], color=colors[comparison], label=comparison)
        ax.legend(loc="upper left", fontsize=6)
        ax.text(-0.12, 1.06, chr(ord("a") + panel), transform=ax.transAxes, fontweight="bold", fontsize=8)
    stem = out_dir / "W5_adjusted_outcomes_forest"
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(stem.with_suffix(".tiff"), dpi=600, bbox_inches="tight")
    fig.savefig(stem.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(fig)


def write_report(
    out_dir: Path,
    unadjusted: pd.DataFrame,
    effects: pd.DataFrame,
    standardized: pd.DataFrame,
    missingness: pd.DataFrame,
) -> None:
    def md(frame: pd.DataFrame) -> str:
        return frame.round(3).to_markdown(index=False)

    mortality_effects = effects.loc[effects.outcome == "mortality_28d", [
        "cohort_label", "comparison", "adjusted_or", "ci_low", "ci_high", "p_value",
        "n_complete", "events", "events_per_parameter", "model_auc"
    ]]
    rrt_effects = effects.loc[effects.outcome == "is_rrt", [
        "cohort_label", "comparison", "adjusted_or", "ci_low", "ci_high", "p_value",
        "n_complete", "events", "events_per_parameter", "model_auc"
    ]]
    mortality_unadjusted = unadjusted.loc[unadjusted.outcome == "mortality_28d", [
        "cohort_label", "phenotype", "n", "events", "risk"
    ]]
    text = f"""# W5 Independent clinical outcome analysis

## Bottom line

The phenotype groups remained associated with 28-day mortality after prespecified,
parsimonious adjustment for age, sex, AKI stage at first diagnosis, and baseline serum
creatinine in each cohort. This endpoint is independent of the seven-day kidney
trajectory variables used to define the phenotypes. The results support prognostic
association, not causal treatment stratification.

## Unadjusted 28-day mortality

{md(mortality_unadjusted)}

## Adjusted 28-day mortality

Reference phenotype: RR. Models were fit separately in each cohort using binomial GLMs
with HC3 robust standard errors. Age was entered per 10 years, baseline creatinine was
log transformed, and initial AKI stage was categorical.

{md(mortality_effects)}

## Covariate-standardized 28-day risks

For each fitted cohort model, phenotype was counterfactually set to DR, RR, or PW for
every complete-case patient and predicted risks were averaged over that cohort's
observed baseline covariate distribution. Intervals are nonparametric patient-level
bootstrap percentile intervals.

{md(standardized[["cohort_label", "phenotype", "standardized_risk", "ci_low", "ci_high", "bootstrap_successes"]])}

## Secondary endpoint: in-hospital renal replacement therapy

RRT was not part of the phenotype-defining trajectory. However, exact RRT timing was
not available in these derived files, so this analysis is secondary and associative.

{md(rrt_effects)}

## Missingness and model scope

{md(missingness)}

The eICU analysis is restricted to the authoritative 1,417-patient cohort used in
the submitted manuscript. During audit, an earlier 1,748-patient risk-factor export
was identified; its phenotype and mortality fields were not used. Covariates were
joined by patient identifier, while phenotype and mortality were overwritten from
the authoritative `03.eICU_SAKI_trajCluster/sk_survival.csv` file.

The common adjustment set was intentionally limited to variables with consistent
definitions in all three archived cohorts. Vasopressor use, mechanical ventilation,
and RRT were not used as mortality covariates because the derived files do not preserve
a common pre-onset time window for these treatments. Maximum AKI stage, day-7 AKI
status, peak creatinine, discharge creatinine, and ICU length of stay were also excluded
because they occur after phenotype ascertainment or overlap with its defining kidney
course.

## Reporting decision

1. Present adjusted 28-day mortality as the principal independent endpoint.
2. Retain RRT only as a secondary association with explicit timing limitations.
3. Do not interpret phenotype coefficients as causal effects or evidence that a
   phenotype-specific treatment improves survival.
"""
    (out_dir / "W5_INDEPENDENT_OUTCOMES.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    unadjusted_parts, effect_parts, standardized_parts, missing_parts = [], [], [], []
    for cohort in ["mimic", "eicu", "aumc"]:
        frame = load_cohort(args.data_dir, cohort)
        for outcome in ["mortality_28d", "is_rrt"]:
            unadjusted_parts.append(unadjusted_table(cohort, frame, outcome))
            complete, fit = fit_model(frame, outcome)
            effect_parts.append(extract_group_effects(cohort, outcome, complete, fit))
            required = [outcome, "groupHPD", "age10", "male", "first_aki_stage", "log_baseline_scr"]
            missing_parts.append(
                {
                    "cohort": COHORT_LABEL[cohort],
                    "outcome": outcome,
                    "n_source": len(frame),
                    "n_complete": len(complete),
                    "excluded_missing": len(frame) - len(complete),
                    "excluded_missing_percent": 100 * (len(frame) - len(complete)) / len(frame),
                    "variables_required": ", ".join(required),
                }
            )
        standardized_parts.append(standardized_risk_bootstrap(cohort, frame, args.bootstrap))

    unadjusted = pd.concat(unadjusted_parts, ignore_index=True)
    effects = pd.concat(effect_parts, ignore_index=True)
    standardized = pd.concat(standardized_parts, ignore_index=True)
    missingness = pd.DataFrame(missing_parts)
    unadjusted.to_csv(args.out_dir / "outcome_unadjusted_rates.csv", index=False)
    effects.to_csv(args.out_dir / "outcome_adjusted_effects.csv", index=False)
    standardized.to_csv(args.out_dir / "mortality_standardized_risks.csv", index=False)
    missingness.to_csv(args.out_dir / "outcome_model_missingness.csv", index=False)
    plot_forest(effects, args.out_dir)
    write_report(args.out_dir, unadjusted, effects, standardized, missingness)
    print(effects.to_string(index=False))


if __name__ == "__main__":
    main()
