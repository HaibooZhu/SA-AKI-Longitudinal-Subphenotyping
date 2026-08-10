#!/usr/bin/env python3
"""Adjusted independent clinical outcome analyses across the three cohorts.

All inputs are read from the frozen snapshot.  Only aggregate results and figures
are written; no patient identifier is exported.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.metrics import roc_auc_score


STYLE_DIR = Path(__file__).resolve().parents[1] / "07_tables_figures"
sys.path.insert(0, str(STYLE_DIR))
from publication_figure_style import (  # noqa: E402
    DOUBLE_COLUMN_IN,
    GRID_LIGHT,
    LINE_AUX,
    PHENOTYPE_COLORS,
    add_panel_label,
    apply_publication_style,
    export_figure,
)


SEED = 20260805
PHENOTYPE = {1: "DR", 2: "RR", 3: "PW"}
COHORT_LABEL = {"mimic": "MIMIC-IV", "eicu": "eICU-CRD", "aumc": "AUMC"}
EXPECTED_N = {"mimic": 4713, "eicu": 1417, "aumc": 2183}
COHORT_FILES = {
    "mimic": ("01.MIMICIV_SAKI_trajCluster", "df_mixAK_fea4_C3.csv"),
    "eicu": ("03.eICU_SAKI_trajCluster", "df_mixAK_fea4_C3_eicu.csv"),
    "aumc": ("02.AUMCdb_SAKI_trajCluster", "df_mixAK_fea3_C3_aumc.csv"),
}
COVARIATE_SOURCES = {
    "mimic": {
        "demographics": "00.data_mimic/feature_data/df_mimic_basicinfo.csv",
        "baseline": "00.data_mimic/disease_definition/AKI/df_base_crea.csv",
        "stage": "00.data_mimic/disease_definition/AKI/sk_first_and_max_stage.csv",
        "rrt": "00.data_mimic/treatment/lifesupport.csv",
        "baseline_id": "stay_id",
        "baseline_column": "baseline_Scr",
        "baseline_multiplier": 1.0,
        "rrt_missing_means_zero": False,
    },
    "eicu": {
        "demographics": "00.data_eicu/feature_data/df_eicu_basicinfo.csv",
        "baseline": "00.data_eicu/disease_definition/AKI/df_base_crea.csv",
        "stage": "00.data_eicu/disease_definition/AKI/eicu_sk_first_and_max_stage.csv",
        "rrt": "00.data_eicu/treatment/eicu_lifesupport.csv",
        "baseline_id": "stay_id",
        "baseline_column": "baseline_creatinine",
        "baseline_multiplier": 1.0,
        "rrt_missing_means_zero": True,
    },
    "aumc": {
        "demographics": "00.data_aumc/feature_data/df_aumc_basicinfo.csv",
        "baseline": "00.data_aumc/disease_definition/AKI/baseline_creatinine.csv",
        "stage": "00.data_aumc/disease_definition/AKI/aumc_first_and_max_stage.csv",
        "rrt": "00.data_aumc/treatment/aumcdb_lifesupport.csv",
        "baseline_id": "admissionid",
        "baseline_column": "baseline_creatinine",
        "baseline_multiplier": 0.01131,
        "rrt_missing_means_zero": False,
    },
}
SOFA_COMPONENTS = (
    "respiration_sofa",
    "coagulation_sofa",
    "liver_sofa",
    "cardiovascular_sofa",
    "cns_sofa",
    "renal_sofa",
)
NONRENAL_SOFA_COMPONENTS = SOFA_COMPONENTS[:-1]
MODEL_FORMULAS = {
    "minimal_baseline": (
        "{outcome} ~ C(groupHPD, Treatment(reference=2)) + age10 + male + "
        "C(first_aki_stage) + log_baseline_scr"
    ),
    "primary_onset_nonrenal_sofa": (
        "{outcome} ~ C(groupHPD, Treatment(reference=2)) + age10 + male + "
        "C(first_aki_stage) + log_baseline_scr + onset_nonrenal_sofa"
    ),
}
MODEL_COLUMNS = {
    "minimal_baseline": (
        "groupHPD",
        "age10",
        "male",
        "first_aki_stage",
        "log_baseline_scr",
    ),
    "primary_onset_nonrenal_sofa": (
        "groupHPD",
        "age10",
        "male",
        "first_aki_stage",
        "log_baseline_scr",
        "onset_nonrenal_sofa",
    ),
}
POPULATION_LABEL = {
    "overall": "Full authoritative cohort",
    "day7_all_survivors": "Alive at day 7",
    "day7_full_trajectory": "Alive at day 7 and observed through time 28",
}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--snapshot-root",
        type=Path,
        default=repo / "00_frozen_inputs/data_snapshot/remote_project_snapshot",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=repo / "02_revision_outputs/reports/W5_independent_outcomes",
    )
    parser.add_argument("--bootstrap", type=int, default=300)
    return parser.parse_args()


def mismatch_count(left: pd.Series, right: pd.Series) -> int:
    """Count disagreements only where both values are observed."""
    comparable = left.notna() & right.notna()
    return int((left.loc[comparable] != right.loc[comparable]).sum())


def collapse_unique_source(
    frame: pd.DataFrame, columns: list[str], source_name: str
) -> pd.DataFrame:
    """Collapse exact duplicate source rows and reject conflicting patient values."""
    selected = frame[columns].dropna(subset=["stay_id"]).copy()
    for column in columns[1:]:
        conflicts = selected.groupby("stay_id")[column].nunique(dropna=True).gt(1)
        if conflicts.any():
            raise ValueError(
                f"{source_name}: conflicting {column} values for "
                f"{int(conflicts.sum())} patients"
            )
    return selected.drop_duplicates("stay_id", keep="first")


def load_upstream_covariates(
    snapshot_root: Path,
    cohort: str,
    authoritative_ids: pd.Series,
    legacy_risk: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Rebuild model covariates from their frozen upstream files.

    The legacy risk-factor exports are compared for audit only. Their values are
    never selected for the revision models.
    """
    spec = COVARIATE_SOURCES[cohort]
    demographics_path = snapshot_root / spec["demographics"]
    baseline_path = snapshot_root / spec["baseline"]
    stage_path = snapshot_root / spec["stage"]
    rrt_path = snapshot_root / spec["rrt"]

    demographics = collapse_unique_source(
        pd.read_csv(demographics_path),
        ["stay_id", "gender", "age"],
        spec["demographics"],
    )
    baseline_raw = pd.read_csv(baseline_path).rename(
        columns={spec["baseline_id"]: "stay_id"}
    )
    baseline = collapse_unique_source(
        baseline_raw,
        ["stay_id", spec["baseline_column"]],
        spec["baseline"],
    ).rename(columns={spec["baseline_column"]: "baseline_Scr"})
    baseline["baseline_Scr"] = (
        pd.to_numeric(baseline["baseline_Scr"], errors="coerce")
        * float(spec["baseline_multiplier"])
    )
    stage = collapse_unique_source(
        pd.read_csv(stage_path),
        ["stay_id", "first_aki_stage"],
        spec["stage"],
    )
    rrt = collapse_unique_source(
        pd.read_csv(rrt_path),
        ["stay_id", "is_rrt"],
        spec["rrt"],
    )

    covariates = pd.DataFrame({"stay_id": authoritative_ids}).drop_duplicates()
    for source in [demographics, baseline, stage, rrt]:
        covariates = covariates.merge(
            source, on="stay_id", how="left", validate="one_to_one"
        )
    if spec["rrt_missing_means_zero"]:
        # This reproduces the explicit eICU notebook rule after the left join.
        covariates["is_rrt"] = covariates["is_rrt"].fillna(0)

    legacy = legacy_risk[
        ["stay_id", "gender", "age", "baseline_Scr", "first_aki_stage", "is_rrt"]
    ].copy()
    audit = covariates.merge(
        legacy,
        on="stay_id",
        how="left",
        suffixes=("_upstream", "_legacy"),
        validate="one_to_one",
    )
    source_by_variable = {
        "gender": spec["demographics"],
        "age": spec["demographics"],
        "baseline_Scr": spec["baseline"],
        "first_aki_stage": spec["stage"],
        "is_rrt": spec["rrt"],
    }
    verification_rows = []
    for variable, source_path in source_by_variable.items():
        upstream = audit[f"{variable}_upstream"]
        legacy_values = audit[f"{variable}_legacy"]
        comparable = upstream.notna() & legacy_values.notna()
        if variable == "gender":
            mismatch = int(
                upstream.loc[comparable]
                .astype(str)
                .ne(legacy_values.loc[comparable].astype(str))
                .sum()
            )
        else:
            left = pd.to_numeric(upstream.loc[comparable], errors="coerce")
            right = pd.to_numeric(legacy_values.loc[comparable], errors="coerce")
            mismatch = int(
                (~np.isclose(left, right, rtol=1e-8, atol=1e-8, equal_nan=True)).sum()
            )
        verification_rows.append(
            {
                "cohort": COHORT_LABEL[cohort],
                "covariate": variable,
                "upstream_source": source_path,
                "authoritative_n": len(covariates),
                "upstream_nonmissing_n": int(upstream.notna().sum()),
                "legacy_nonmissing_n": int(legacy_values.notna().sum()),
                "comparable_n": int(comparable.sum()),
                "legacy_vs_upstream_mismatch_n": mismatch,
                "revision_value_source": "frozen upstream file",
                "missing_value_rule": (
                    "left-join absence explicitly set to 0 in archived eICU notebook"
                    if variable == "is_rrt" and spec["rrt_missing_means_zero"]
                    else "preserved as missing"
                ),
            }
        )
    return covariates, pd.DataFrame(verification_rows)


def load_authoritative(snapshot_root: Path, cohort: str) -> pd.DataFrame:
    cohort_dir, _ = COHORT_FILES[cohort]
    path = snapshot_root / cohort_dir / "sk_survival.csv"
    authoritative = pd.read_csv(path)[
        ["stay_id", "groupHPD", "mortality_28d", "mortality_7d"]
    ].copy()
    if authoritative.stay_id.duplicated().any():
        raise ValueError(f"Duplicate authoritative stay_id values in {cohort}")
    if len(authoritative) != EXPECTED_N[cohort]:
        raise ValueError(
            f"Expected {EXPECTED_N[cohort]:,} authoritative {cohort} patients, "
            f"found {len(authoritative):,}"
        )
    for column in ["groupHPD", "mortality_28d", "mortality_7d"]:
        authoritative[column] = pd.to_numeric(authoritative[column], errors="coerce")
    if not authoritative.groupHPD.dropna().isin(PHENOTYPE).all():
        raise ValueError(f"Unexpected authoritative phenotype label in {cohort}")
    return authoritative


def load_onset_nonrenal_sofa(snapshot_root: Path, cohort: str) -> pd.DataFrame:
    path = (
        snapshot_root
        / "04.other_feature_in_three_dataset/03.sofa_feature"
        / f"{cohort}_sofa_clean.csv"
    )
    sofa = pd.read_csv(path, usecols=["stay_id", "time", *SOFA_COMPONENTS])
    sofa = sofa.loc[pd.to_numeric(sofa.time, errors="coerce").eq(1)].copy()
    if sofa.stay_id.duplicated().any():
        raise ValueError(f"Duplicate onset-window SOFA rows in {cohort}")
    for column in SOFA_COMPONENTS:
        sofa[column] = pd.to_numeric(sofa[column], errors="coerce")
        if sofa[column].dropna().lt(0).any() or sofa[column].dropna().gt(4).any():
            raise ValueError(f"Out-of-range {column} values in {cohort}")
    # Deliberately exclude the renal component because the phenotype is defined
    # from creatinine and urine-output trajectories.
    sofa["onset_nonrenal_sofa"] = sofa[list(NONRENAL_SOFA_COMPONENTS)].sum(
        axis=1, min_count=len(NONRENAL_SOFA_COMPONENTS)
    )
    return sofa[["stay_id", "onset_nonrenal_sofa"]]


def load_trajectory_summary(snapshot_root: Path, cohort: str) -> pd.DataFrame:
    cohort_dir, trajectory_file = COHORT_FILES[cohort]
    trajectory = pd.read_csv(
        snapshot_root / cohort_dir / trajectory_file, usecols=["stay_id", "time"]
    )
    trajectory["time"] = pd.to_numeric(trajectory.time, errors="coerce")
    if trajectory[["stay_id", "time"]].duplicated().any():
        raise ValueError(f"Duplicate patient-time rows in {cohort} trajectory")
    summary = (
        trajectory.groupby("stay_id", as_index=False)
        .agg(
            trajectory_min_time=("time", "min"),
            trajectory_max_time=("time", "max"),
            trajectory_rows=("time", "size"),
        )
    )
    summary["observed_through_time28"] = summary.trajectory_max_time.ge(28)
    return summary


def load_cohort(
    snapshot_root: Path, cohort: str
) -> tuple[pd.DataFrame, dict, pd.DataFrame]:
    risk_dir = (
        snapshot_root
        / "04.other_feature_in_three_dataset/08.subphenotype_association_analysis"
    )
    files = [risk_dir / f"df_{cohort}_c{group}_riskfactor.csv" for group in [1, 2, 3]]
    risk = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    if risk.stay_id.duplicated().any():
        raise ValueError(f"Duplicate stay_id values in {cohort} risk-factor export")

    authoritative = load_authoritative(snapshot_root, cohort)
    overlap = authoritative.merge(
        risk[["stay_id", "groupHPD", "mortality_28d"]],
        on="stay_id",
        how="left",
        suffixes=("_authoritative", "_legacy"),
        validate="one_to_one",
        indicator=True,
    )
    missing_covariate_rows = int(overlap._merge.ne("both").sum())
    if missing_covariate_rows:
        raise ValueError(
            f"Risk-factor export misses {missing_covariate_rows} authoritative {cohort} patients"
        )
    lineage = {
        "cohort": COHORT_LABEL[cohort],
        "authoritative_n": len(authoritative),
        "legacy_risk_n": len(risk),
        "legacy_extra_n": int((~risk.stay_id.isin(authoritative.stay_id)).sum()),
        "legacy_group_mismatch_n": mismatch_count(
            overlap.groupHPD_authoritative, overlap.groupHPD_legacy
        ),
        "legacy_mortality_mismatch_n": mismatch_count(
            overlap.mortality_28d_authoritative, overlap.mortality_28d_legacy
        ),
        "authoritative_mortality28_missing_n": int(
            authoritative.mortality_28d.isna().sum()
        ),
        "authoritative_mortality7_missing_n": int(
            authoritative.mortality_7d.isna().sum()
        ),
    }

    covariates, covariate_verification = load_upstream_covariates(
        snapshot_root, cohort, authoritative["stay_id"], risk
    )
    frame = authoritative.merge(
        covariates, on="stay_id", how="left", validate="one_to_one"
    )
    sofa = load_onset_nonrenal_sofa(snapshot_root, cohort)
    trajectory = load_trajectory_summary(snapshot_root, cohort)
    frame = frame.merge(sofa, on="stay_id", how="left", validate="one_to_one")
    frame = frame.merge(trajectory, on="stay_id", how="left", validate="one_to_one")
    if frame.onset_nonrenal_sofa.isna().any():
        raise ValueError(f"Missing onset nonrenal SOFA for authoritative {cohort} patients")
    if frame.trajectory_max_time.isna().any():
        raise ValueError(f"Missing trajectory summary for authoritative {cohort} patients")

    frame["groupHPD"] = pd.to_numeric(frame.groupHPD, errors="coerce").astype("Int64")
    frame["age10"] = pd.to_numeric(frame.age, errors="coerce") / 10
    frame["male"] = frame.gender.map({"M": 1.0, "F": 0.0})
    baseline = pd.to_numeric(frame.baseline_Scr, errors="coerce")
    frame["log_baseline_scr"] = np.where(baseline > 0, np.log(baseline), np.nan)
    for column in ["first_aki_stage", "mortality_28d", "mortality_7d", "is_rrt"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["mortality_day8_28"] = frame.mortality_28d.where(frame.mortality_7d.eq(0))
    frame["landmark_day7_survivor"] = frame.mortality_7d.eq(0)
    frame["landmark_full_trajectory"] = (
        frame.landmark_day7_survivor & frame.observed_through_time28
    )
    return frame, lineage, covariate_verification


def select_population(frame: pd.DataFrame, population: str) -> pd.DataFrame:
    if population == "overall":
        return frame.copy()
    if population == "day7_all_survivors":
        return frame.loc[frame.landmark_day7_survivor].copy()
    if population == "day7_full_trajectory":
        return frame.loc[frame.landmark_full_trajectory].copy()
    raise ValueError(f"Unknown population: {population}")


def fit_model(
    frame: pd.DataFrame,
    outcome: str,
    model_variant: str,
    population: str = "overall",
):
    selected = select_population(frame, population)
    columns = [outcome, *MODEL_COLUMNS[model_variant]]
    complete = selected[columns].dropna().copy()
    complete["groupHPD"] = complete.groupHPD.astype(int)
    complete["first_aki_stage"] = complete.first_aki_stage.astype(int)
    if complete[outcome].nunique() != 2:
        raise ValueError(
            f"Outcome {outcome} does not contain both classes in {population}"
        )
    formula = MODEL_FORMULAS[model_variant].format(outcome=outcome)
    fit = smf.glm(formula, data=complete, family=sm.families.Binomial()).fit(
        cov_type="HC3"
    )
    return selected, complete, fit


def extract_group_effects(
    cohort: str,
    outcome: str,
    model_variant: str,
    population: str,
    selected: pd.DataFrame,
    complete: pd.DataFrame,
    fit,
) -> pd.DataFrame:
    rows = []
    event_count = int(complete[outcome].sum())
    n_parameters = len(fit.params)
    group_events = complete.groupby("groupHPD")[outcome].sum()
    for group in [1, 3]:
        term = f"C(groupHPD, Treatment(reference=2))[T.{group}]"
        beta = fit.params[term]
        low, high = fit.conf_int().loc[term]
        rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABEL[cohort],
                "outcome": outcome,
                "analysis_population": population,
                "population_label": POPULATION_LABEL[population],
                "model_variant": model_variant,
                "comparison": f"{PHENOTYPE[group]} vs RR",
                "group": group,
                "reference_group": 2,
                "adjusted_or": np.exp(beta),
                "ci_low": np.exp(low),
                "ci_high": np.exp(high),
                "p_value": fit.pvalues[term],
                "n_population": len(selected),
                "n_complete": len(complete),
                "events": event_count,
                "events_in_compared_group": int(group_events.get(group, 0)),
                "events_in_reference_group": int(group_events.get(2, 0)),
                "n_parameters": n_parameters,
                "events_per_parameter": float(event_count / n_parameters),
                "model_auc": roc_auc_score(complete[outcome], fit.predict(complete)),
                "converged": bool(fit.converged),
                "diagnostic_flag": (
                    "PASS" if event_count / n_parameters >= 10 else "CAUTION_LOW_EVENTS"
                ),
            }
        )
    return pd.DataFrame(rows)


def standardized_risk_bootstrap(
    cohort: str,
    frame: pd.DataFrame,
    outcome: str,
    population: str,
    model_variant: str,
    n_boot: int,
) -> pd.DataFrame:
    _, complete, fit = fit_model(frame, outcome, model_variant, population)
    point = {}
    for group in [1, 2, 3]:
        counterfactual = complete.copy()
        counterfactual["groupHPD"] = group
        point[group] = float(np.mean(fit.predict(counterfactual)))

    population_seed = {"overall": 0, "day7_full_trajectory": 100}[population]
    cohort_seed = {"mimic": 1, "eicu": 2, "aumc": 3}[cohort]
    rng = np.random.default_rng(SEED + population_seed + cohort_seed)
    draws = {group: [] for group in [1, 2, 3]}
    failures = 0
    formula = MODEL_FORMULAS[model_variant].format(outcome=outcome)
    for _ in range(n_boot):
        boot = complete.iloc[rng.integers(0, len(complete), len(complete))].copy()
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                boot_fit = smf.glm(
                    formula, data=boot, family=sm.families.Binomial()
                ).fit()
            for group in [1, 2, 3]:
                counterfactual = complete.copy()
                counterfactual["groupHPD"] = group
                prediction = float(np.mean(boot_fit.predict(counterfactual)))
                if not np.isfinite(prediction):
                    raise ValueError("Non-finite bootstrap prediction")
                draws[group].append(prediction)
        except Exception:
            failures += 1

    rows = []
    for group in [1, 2, 3]:
        if not draws[group]:
            raise RuntimeError(f"All bootstrap fits failed for {cohort} {population}")
        rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABEL[cohort],
                "outcome": outcome,
                "analysis_population": population,
                "model_variant": model_variant,
                "phenotype": PHENOTYPE[group],
                "group": group,
                "standardized_risk": point[group],
                "ci_low": np.quantile(draws[group], 0.025),
                "ci_high": np.quantile(draws[group], 0.975),
                "bootstrap_successes": len(draws[group]),
                "bootstrap_failures": failures,
            }
        )
    return pd.DataFrame(rows)


def unadjusted_table(
    cohort: str, frame: pd.DataFrame, outcome: str, population: str
) -> pd.DataFrame:
    selected = select_population(frame, population)
    work = selected[["groupHPD", outcome]].dropna().copy()
    work["groupHPD"] = work.groupHPD.astype(int)
    result = (
        work.groupby("groupHPD")[outcome]
        .agg(n="size", events="sum", risk="mean")
        .reset_index()
    )
    result.insert(0, "cohort", cohort)
    result.insert(1, "cohort_label", COHORT_LABEL[cohort])
    result.insert(2, "outcome", outcome)
    result.insert(3, "analysis_population", population)
    result["phenotype"] = result.groupHPD.map(PHENOTYPE)
    return result


def population_table(cohort: str, frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for population in ["overall", "day7_all_survivors", "day7_full_trajectory"]:
        selected = select_population(frame, population)
        outcome = "mortality_28d" if population == "overall" else "mortality_day8_28"
        observed = selected[outcome].notna()
        rows.append(
            {
                "cohort": cohort,
                "cohort_label": COHORT_LABEL[cohort],
                "analysis_population": population,
                "population_label": POPULATION_LABEL[population],
                "n_selected": len(selected),
                "n_with_outcome": int(observed.sum()),
                "events": int(selected.loc[observed, outcome].sum()),
                "event_risk": float(selected.loc[observed, outcome].mean()),
                "retained_from_authoritative": len(selected) / len(frame),
            }
        )
    return pd.DataFrame(rows)


def plot_forest(effects: pd.DataFrame, out_dir: Path) -> None:
    apply_publication_style()
    panels = [
        ("mortality_28d", "overall", "28-day mortality"),
        (
            "mortality_day8_28",
            "day7_full_trajectory",
            "Day 8–28 mortality: strict day-7 landmark",
        ),
    ]
    colors = {
        "DR vs RR": PHENOTYPE_COLORS["DR"],
        "PW vs RR": PHENOTYPE_COLORS["PW"],
    }
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(DOUBLE_COLUMN_IN, 3.35),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    y_positions = {
        ("MIMIC-IV", "DR vs RR"): 5,
        ("MIMIC-IV", "PW vs RR"): 4,
        ("eICU-CRD", "DR vs RR"): 3,
        ("eICU-CRD", "PW vs RR"): 2,
        ("AUMC", "DR vs RR"): 1,
        ("AUMC", "PW vs RR"): 0,
    }
    for panel, (outcome, population, title) in enumerate(panels):
        ax = axes[panel]
        subset = effects.loc[
            effects.outcome.eq(outcome)
            & effects.analysis_population.eq(population)
            & effects.model_variant.eq("primary_onset_nonrenal_sofa")
        ]
        for _, row in subset.iterrows():
            y = y_positions[(row.cohort_label, row.comparison)]
            ax.errorbar(
                row.adjusted_or,
                y,
                xerr=[[row.adjusted_or - row.ci_low], [row.ci_high - row.adjusted_or]],
                fmt="o",
                markersize=4.2,
                capsize=2.0,
                linewidth=1.1,
                color=colors[row.comparison],
            )
        for y0, y1 in [(3.55, 5.45), (-0.45, 1.45)]:
            ax.axhspan(y0, y1, color="#F6F6F6", zorder=0)
        ax.axvline(1, color="#777777", linestyle="--", linewidth=LINE_AUX)
        ax.set_xscale("log")
        ax.set_xlim(0.88, 16)
        ax.set_ylim(-0.55, 5.55)
        ax.set_yticks(
            [5, 4, 3, 2, 1, 0],
            [
                "MIMIC-IV · DR",
                "MIMIC-IV · PW",
                "eICU-CRD · DR",
                "eICU-CRD · PW",
                "AUMC · DR",
                "AUMC · PW",
            ],
        )
        ax.set_xlabel("Adjusted odds ratio (95% CI)")
        ax.set_title(title, fontweight="bold", pad=5)
        ax.grid(axis="x", color=GRID_LIGHT, linewidth=0.5)
        add_panel_label(ax, chr(ord("a") + panel))
    stem = out_dir / "W5_adjusted_outcomes_forest"
    export_figure(fig, stem)
    plt.close(fig)


def write_report(
    out_dir: Path,
    lineage: pd.DataFrame,
    covariate_verification: pd.DataFrame,
    populations: pd.DataFrame,
    unadjusted: pd.DataFrame,
    effects: pd.DataFrame,
    standardized: pd.DataFrame,
    missingness: pd.DataFrame,
) -> None:
    def md(frame: pd.DataFrame) -> str:
        return frame.round(3).to_markdown(index=False)

    display_effect_columns = [
        "cohort_label",
        "comparison",
        "adjusted_or",
        "ci_low",
        "ci_high",
        "p_value",
        "n_complete",
        "events",
        "events_per_parameter",
        "model_auc",
        "diagnostic_flag",
    ]
    overall_primary = effects.loc[
        effects.outcome.eq("mortality_28d")
        & effects.analysis_population.eq("overall")
        & effects.model_variant.eq("primary_onset_nonrenal_sofa"),
        display_effect_columns,
    ]
    overall_minimal = effects.loc[
        effects.outcome.eq("mortality_28d")
        & effects.analysis_population.eq("overall")
        & effects.model_variant.eq("minimal_baseline"),
        display_effect_columns,
    ]
    landmark_strict = effects.loc[
        effects.outcome.eq("mortality_day8_28")
        & effects.analysis_population.eq("day7_full_trajectory"),
        display_effect_columns,
    ]
    landmark_all = effects.loc[
        effects.outcome.eq("mortality_day8_28")
        & effects.analysis_population.eq("day7_all_survivors"),
        display_effect_columns,
    ]
    rrt_effects = effects.loc[
        effects.outcome.eq("is_rrt"), display_effect_columns
    ]
    mortality_unadjusted = unadjusted.loc[
        unadjusted.outcome.eq("mortality_28d"),
        ["cohort_label", "phenotype", "n", "events", "risk"],
    ]
    landmark_unadjusted = unadjusted.loc[
        unadjusted.outcome.eq("mortality_day8_28"),
        ["cohort_label", "analysis_population", "phenotype", "n", "events", "risk"],
    ]
    standard_columns = [
        "cohort_label",
        "analysis_population",
        "phenotype",
        "standardized_risk",
        "ci_low",
        "ci_high",
        "bootstrap_successes",
        "bootstrap_failures",
    ]
    text = f"""# W5 Independent clinical outcome analysis

## Bottom line

This revision analysis uses authoritative phenotype and survival files for every
cohort and adjusts mortality associations for age, sex, AKI stage at first diagnosis,
baseline serum creatinine, and onset-window nonrenal SOFA. The endpoint is independent
of the kidney variables used to define the seven-day phenotypes. Results support
prognostic association only; they do not establish causal treatment response.

The strict day-7 landmark analysis is a sensitivity analysis. Requiring a complete
trajectory through time 28 preferentially retains longer-observed ICU stays and can
introduce selection/collider bias. Therefore, the all-day-7-survivor analysis is shown
alongside it, and neither replaces the full-cohort mortality analysis.

## Authoritative cohort and lineage checks

Phenotype, 28-day mortality, and 7-day mortality were overwritten from each cohort's
frozen `sk_survival.csv`. Age, sex, baseline serum creatinine, first AKI stage,
and RRT were rebuilt directly from their frozen upstream source files. The legacy
risk-factor exports were used only as comparison targets and did not supply model
values.

{md(lineage)}

Aggregate upstream-versus-legacy covariate verification:

{md(covariate_verification)}

The 1,748-row legacy eICU risk-factor export is not an analysis population. It contains
331 patients outside the authoritative 1,417-patient cohort; its legacy phenotype and
mortality columns were ignored.

## Outcome populations

{md(populations)}

`day7_all_survivors` includes all patients documented alive at day 7 with an observed
28-day outcome. `day7_full_trajectory` additionally requires a kidney trajectory row
at time 28; this is a restrictive sensitivity population, not a less-biased primary
cohort.

## Covariate definition

The primary model is:

`outcome ~ phenotype + age/10 + sex + first AKI stage + log(baseline SCr) + onset nonrenal SOFA`

Onset nonrenal SOFA is the sum of respiratory, coagulation, liver, cardiovascular, and
CNS SOFA components at aligned time 1 (the first six-hour window after SA-AKI onset).
The renal component is deliberately excluded to avoid readjusting for creatinine or
urine output, which define the phenotypes. These archived SOFA components were derived
from rolling 24-hour scores, aligned to six-hour windows, aggregated by the within-window
maximum, and completed within patient by forward/backward filling; residual missing
values were set to zero in AUMC and eICU. This provenance limits causal interpretation.

## Unadjusted 28-day mortality

{md(mortality_unadjusted)}

## Primary adjusted 28-day mortality

Reference phenotype: RR. Cohort-specific binomial GLMs use HC3 robust standard errors.

{md(overall_primary)}

## Minimal-model sensitivity

This model omits nonrenal SOFA and otherwise uses the same baseline covariates.

{md(overall_minimal)}

## Day-7 landmark analyses: death during days 8–28

Unadjusted event rates:

{md(landmark_unadjusted)}

Strict full-trajectory sensitivity:

{md(landmark_strict)}

All documented day-7 survivors:

{md(landmark_all)}

## Covariate-standardized mortality risks

Phenotype was counterfactually set to DR, RR, or PW for every complete-case patient;
predictions were averaged over that cohort's observed covariate distribution. Intervals
are patient-level nonparametric bootstrap percentile intervals.

{md(standardized[standard_columns])}

## Secondary endpoint: in-hospital renal replacement therapy

RRT was not part of the phenotype-defining trajectory, but exact treatment timing is
not available in these derived files. It is therefore secondary and associative.

{md(rrt_effects)}

## Missingness and diagnostics

{md(missingness)}

`CAUTION_LOW_EVENTS` means fewer than 10 outcome events per fitted coefficient. It is a
warning about precision and possible overfitting, not an automatic model failure.

## Reporting decision

1. Use the nonrenal-SOFA-adjusted full-cohort 28-day mortality model as the principal
   independent clinical endpoint.
2. Report both day-7 landmark definitions and explicitly label the full-trajectory
   version as a selection-sensitive analysis.
3. Retain RRT only as a secondary association with a treatment-timing limitation.
4. Do not interpret phenotype coefficients as causal effects or evidence that a
   phenotype-specific treatment improves survival.
"""
    (out_dir / "W5_INDEPENDENT_OUTCOMES.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    lineage_rows = []
    covariate_verification_parts = []
    population_parts = []
    unadjusted_parts = []
    effect_parts = []
    standardized_parts = []
    missing_parts = []

    for cohort in ["mimic", "eicu", "aumc"]:
        frame, lineage, covariate_verification = load_cohort(args.snapshot_root, cohort)
        lineage_rows.append(lineage)
        covariate_verification_parts.append(covariate_verification)
        population_parts.append(population_table(cohort, frame))

        unadjusted_parts.append(unadjusted_table(cohort, frame, "mortality_28d", "overall"))
        for population in ["day7_all_survivors", "day7_full_trajectory"]:
            unadjusted_parts.append(
                unadjusted_table(cohort, frame, "mortality_day8_28", population)
            )

        specifications = [
            ("mortality_28d", "minimal_baseline", "overall"),
            ("mortality_28d", "primary_onset_nonrenal_sofa", "overall"),
            (
                "mortality_day8_28",
                "primary_onset_nonrenal_sofa",
                "day7_all_survivors",
            ),
            (
                "mortality_day8_28",
                "primary_onset_nonrenal_sofa",
                "day7_full_trajectory",
            ),
            ("is_rrt", "primary_onset_nonrenal_sofa", "overall"),
        ]
        for outcome, model_variant, population in specifications:
            selected, complete, fit = fit_model(
                frame, outcome, model_variant, population
            )
            effect_parts.append(
                extract_group_effects(
                    cohort,
                    outcome,
                    model_variant,
                    population,
                    selected,
                    complete,
                    fit,
                )
            )
            required = [outcome, *MODEL_COLUMNS[model_variant]]
            missing_parts.append(
                {
                    "cohort": COHORT_LABEL[cohort],
                    "outcome": outcome,
                    "analysis_population": population,
                    "model_variant": model_variant,
                    "n_population": len(selected),
                    "n_complete": len(complete),
                    "excluded_missing": len(selected) - len(complete),
                    "excluded_missing_percent": (
                        100 * (len(selected) - len(complete)) / len(selected)
                        if len(selected)
                        else np.nan
                    ),
                    "variables_required": ", ".join(required),
                }
            )

        for outcome, population in [
            ("mortality_28d", "overall"),
            ("mortality_day8_28", "day7_full_trajectory"),
        ]:
            standardized_parts.append(
                standardized_risk_bootstrap(
                    cohort,
                    frame,
                    outcome,
                    population,
                    "primary_onset_nonrenal_sofa",
                    args.bootstrap,
                )
            )

    lineage = pd.DataFrame(lineage_rows)
    covariate_verification = pd.concat(
        covariate_verification_parts, ignore_index=True
    )
    populations = pd.concat(population_parts, ignore_index=True)
    unadjusted = pd.concat(unadjusted_parts, ignore_index=True)
    effects = pd.concat(effect_parts, ignore_index=True)
    standardized = pd.concat(standardized_parts, ignore_index=True)
    missingness = pd.DataFrame(missing_parts)

    lineage.to_csv(args.out_dir / "outcome_source_lineage.csv", index=False)
    covariate_verification.to_csv(
        args.out_dir / "covariate_source_verification.csv", index=False
    )
    populations.to_csv(args.out_dir / "outcome_populations.csv", index=False)
    unadjusted.to_csv(args.out_dir / "outcome_unadjusted_rates.csv", index=False)
    effects.to_csv(args.out_dir / "outcome_adjusted_effects.csv", index=False)
    standardized.to_csv(args.out_dir / "mortality_standardized_risks.csv", index=False)
    missingness.to_csv(args.out_dir / "outcome_model_missingness.csv", index=False)
    plot_forest(effects, args.out_dir)
    write_report(
        args.out_dir,
        lineage,
        covariate_verification,
        populations,
        unadjusted,
        effects,
        standardized,
        missingness,
    )
    print(
        effects[
            [
                "cohort_label",
                "outcome",
                "analysis_population",
                "model_variant",
                "comparison",
                "adjusted_or",
                "ci_low",
                "ci_high",
                "p_value",
                "n_complete",
                "events",
                "diagnostic_flag",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
