"""W20 -- Informative-selection audit of the documentation-density restrictions.

Two sensitivity analyses requested during revision restrict the analysed
population by how densely a patient is documented:

  * the >=50% urine-output coverage subset (eICU-CRD), and
  * the complete 30-window follow-up subset (all three cohorts).

Both restrictions were reported as "no stable Progressive Worsening component".
This script asks the prior question the reports did not answer: *who* the
restrictions remove.  If documentation density is independent of phenotype and
of survival, each restriction should retain the phenotypes in proportion.  If it
is not independent, the restricted population cannot test whether the phenotype
exists, because the patients who define the phenotype are gone before the model
is fitted.

All inputs are frozen; only aggregate results are written.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
FROZEN = REPO / "00_frozen_inputs/data_snapshot"
REMOTE = FROZEN / "remote_project_snapshot"
LOCAL = FROZEN / "local_derived_data"
COVERAGE = (
    REPO
    / "02_revision_outputs/analysis_inputs/eicu_uo_sensitivity/eicu_uo_patient_coverage.csv"
)
OUT = REPO / "02_revision_outputs/reports/W20_documentation_selection"

LABELS = {1: "DR", 2: "RR", 3: "PW"}
COHORTS = {
    "MIMIC-IV": {
        "survival": REMOTE / "01.MIMICIV_SAKI_trajCluster/sk_survival.csv",
        "matrix": LOCAL / "df_mixAK_fea4_C3_mimic.csv",
        "features": ["bun", "creatinine", "urineoutput", "crea_divide_basecrea"],
    },
    "eICU-CRD": {
        "survival": REMOTE / "03.eICU_SAKI_trajCluster/sk_survival.csv",
        "matrix": LOCAL / "df_mixAK_fea4_C3_eicu.csv",
        "features": ["bun", "creatinine", "urineoutput", "crea_divide_basecrea"],
    },
    "AUMC": {
        "survival": REMOTE / "02.AUMCdb_SAKI_trajCluster/sk_survival.csv",
        "matrix": LOCAL / "df_mixAK_fea3_C3_aumc.csv",
        "features": ["creatinine", "urineoutput", "crea_divide_basecrea"],
    },
}
PLANNED_WINDOWS = 30


def _phenotype(frame: pd.DataFrame) -> pd.Series:
    return frame["groupHPD"].astype(int).map(LABELS)


def _retention_table(
    frame: pd.DataFrame, keep: pd.Series, cohort: str, restriction: str
) -> pd.DataFrame:
    """Retention of each phenotype, and mortality of kept versus removed."""
    rows = []
    for group in ["RR", "DR", "PW", "ALL"]:
        mask = frame["phenotype"].eq(group) if group != "ALL" else pd.Series(
            True, index=frame.index
        )
        sub = frame[mask]
        kept = keep[mask]
        n = int(len(sub))
        n_kept = int(kept.sum())
        mort_all = sub["mortality_28d"]
        rows.append(
            {
                "cohort": cohort,
                "restriction": restriction,
                "phenotype": group,
                "n_total": n,
                "n_retained": n_kept,
                "retention_rate": n_kept / n if n else np.nan,
                "mortality28_total": mort_all.mean(),
                "mortality28_retained": sub.loc[kept.values, "mortality_28d"].mean()
                if n_kept
                else np.nan,
                "mortality28_removed": sub.loc[~kept.values, "mortality_28d"].mean()
                if n - n_kept
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    retention_frames: list[pd.DataFrame] = []
    notes: dict[str, object] = {}

    # ---------- Part A: eICU urine-output documentation coverage ----------
    coverage = pd.read_csv(COVERAGE)
    survival = pd.read_csv(COHORTS["eICU-CRD"]["survival"])
    eicu = coverage.merge(
        survival[["stay_id", "mortality_28d", "mortality_7d"]], on="stay_id", how="left"
    )
    eicu["phenotype"] = _phenotype(eicu)
    if eicu["mortality_28d"].isna().any():
        notes["eicu_mortality28_missing_n"] = int(eicu["mortality_28d"].isna().sum())
        eicu = eicu.dropna(subset=["mortality_28d"])

    coverage_by_phenotype = (
        eicu.groupby("phenotype")["documented_fraction_planned"]
        .agg(
            n="count",
            mean="mean",
            median="median",
            q25=lambda s: s.quantile(0.25),
            q75=lambda s: s.quantile(0.75),
        )
        .reindex(["RR", "DR", "PW"])
        .round(4)
    )
    coverage_by_phenotype["median_documented_windows"] = (
        eicu.groupby("phenotype")["documented_windows"].median().reindex(["RR", "DR", "PW"])
    )
    coverage_by_phenotype.to_csv(OUT / "eicu_coverage_by_phenotype.csv")

    coverage_by_outcome = (
        eicu.groupby(eicu["mortality_28d"].astype(int))["documented_fraction_planned"]
        .agg(n="count", mean="mean", median="median")
        .round(4)
    )
    coverage_by_outcome.index = ["28-day survivor", "28-day decedent"]
    coverage_by_outcome.to_csv(OUT / "eicu_coverage_by_outcome.csv")

    high_coverage = eicu["documented_windows"] >= 15
    retention_frames.append(
        _retention_table(eicu, high_coverage, "eICU-CRD", ">=50% urine-output coverage")
    )

    # Composition shift induced by the coverage filter.
    composition = pd.DataFrame(
        {
            "full_cohort": eicu["phenotype"].value_counts(normalize=True),
            "high_coverage_subset": eicu.loc[high_coverage, "phenotype"].value_counts(
                normalize=True
            ),
        }
    ).reindex(["RR", "DR", "PW"]).round(4)
    composition["n_full"] = eicu["phenotype"].value_counts().reindex(["RR", "DR", "PW"])
    composition["n_high_coverage"] = (
        eicu.loc[high_coverage, "phenotype"].value_counts().reindex(["RR", "DR", "PW"])
    )
    composition.to_csv(OUT / "eicu_high_coverage_composition.csv")

    # Rank correlation between coverage and phenotype severity, and with death.
    severity = eicu["phenotype"].map({"RR": 0, "DR": 1, "PW": 2})
    notes["eicu_spearman_coverage_vs_severity"] = float(
        pd.Series(eicu["documented_fraction_planned"]).corr(severity, method="spearman")
    )
    notes["eicu_spearman_coverage_vs_mortality28"] = float(
        pd.Series(eicu["documented_fraction_planned"]).corr(
            eicu["mortality_28d"].astype(float), method="spearman"
        )
    )

    # ---------- Part B: complete 30-window follow-up ----------
    for cohort, spec in COHORTS.items():
        matrix = pd.read_csv(spec["matrix"])
        required = ["stay_id", "time", *spec["features"]]
        matrix = matrix.dropna(subset=required)[required]
        windows = matrix.groupby("stay_id")["time"].nunique().rename("available_windows")
        surv = pd.read_csv(spec["survival"])[["stay_id", "groupHPD", "mortality_28d"]]
        frame = surv.merge(windows, on="stay_id", how="left")
        frame["available_windows"] = frame["available_windows"].fillna(0)
        frame["phenotype"] = _phenotype(frame)
        frame = frame.dropna(subset=["mortality_28d"])
        keep = frame["available_windows"] >= PLANNED_WINDOWS
        retention_frames.append(
            _retention_table(frame, keep, cohort, "complete 30-window follow-up")
        )

    retention = pd.concat(retention_frames, ignore_index=True)
    numeric = retention.select_dtypes(include="number").columns
    retention[numeric] = retention[numeric].round(4)
    retention.to_csv(OUT / "restriction_retention.csv", index=False)

    with open(OUT / "W20_notes.json", "w") as handle:
        json.dump(notes, handle, indent=2)

    pd.set_option("display.width", 200)
    print("=== eICU urine-output coverage by phenotype ===")
    print(coverage_by_phenotype)
    print("\n=== eICU urine-output coverage by 28-day outcome ===")
    print(coverage_by_outcome)
    print("\n=== eICU >=50% coverage subset composition ===")
    print(composition)
    print("\n=== Retention under each documentation restriction ===")
    print(retention.to_string(index=False))
    print("\n=== notes ===")
    print(json.dumps(notes, indent=2))


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------
# Appended: formal tests and figure for the informative-selection mechanism.
# ---------------------------------------------------------------------------


def selection_tests() -> None:
    """Test whether each restriction removes patients independently of death."""
    from scipy import stats

    retention = pd.read_csv(OUT / "restriction_retention.csv")
    coverage = pd.read_csv(COVERAGE)
    rows = []

    def _one(frame: pd.DataFrame, keep: pd.Series, cohort: str, restriction: str) -> None:
        for group in ["RR", "DR", "PW", "ALL"]:
            mask = (
                frame["phenotype"].eq(group)
                if group != "ALL"
                else pd.Series(True, index=frame.index)
            )
            sub = frame[mask]
            kept = keep[mask].values
            died = sub["mortality_28d"].astype(int).values
            table = np.array(
                [
                    [int(((~kept) & (died == 1)).sum()), int(((~kept) & (died == 0)).sum())],
                    [int((kept & (died == 1)).sum()), int((kept & (died == 0)).sum())],
                ]
            )
            if table.min() < 0 or table.sum(axis=1).min() == 0 or table.sum(axis=0).min() == 0:
                odds, p = np.nan, np.nan
            else:
                odds, p = stats.fisher_exact(table)
            rows.append(
                {
                    "cohort": cohort,
                    "restriction": restriction,
                    "phenotype": group,
                    "removed_deaths": table[0, 0],
                    "removed_survivors": table[0, 1],
                    "retained_deaths": table[1, 0],
                    "retained_survivors": table[1, 1],
                    # fisher_exact on [[removed_deaths, removed_survivors],
                    # [retained_deaths, retained_survivors]] already returns the
                    # odds of death among removed relative to retained patients.
                    "odds_ratio_death_given_removed": odds,
                    "fisher_p": p,
                }
            )

    survival = pd.read_csv(COHORTS["eICU-CRD"]["survival"])
    eicu = coverage.merge(
        survival[["stay_id", "mortality_28d"]], on="stay_id", how="left"
    ).dropna(subset=["mortality_28d"])
    eicu["phenotype"] = _phenotype(eicu)
    _one(eicu, eicu["documented_windows"] >= 15, "eICU-CRD", ">=50% urine-output coverage")

    for cohort, spec in COHORTS.items():
        matrix = pd.read_csv(spec["matrix"])
        required = ["stay_id", "time", *spec["features"]]
        matrix = matrix.dropna(subset=required)[required]
        windows = matrix.groupby("stay_id")["time"].nunique().rename("available_windows")
        surv = pd.read_csv(spec["survival"])[["stay_id", "groupHPD", "mortality_28d"]]
        frame = surv.merge(windows, on="stay_id", how="left")
        frame["available_windows"] = frame["available_windows"].fillna(0)
        frame["phenotype"] = _phenotype(frame)
        frame = frame.dropna(subset=["mortality_28d"]).reset_index(drop=True)
        _one(
            frame,
            frame["available_windows"] >= PLANNED_WINDOWS,
            cohort,
            "complete 30-window follow-up",
        )

    tests = pd.DataFrame(rows)
    tests["fisher_p"] = tests["fisher_p"].map(lambda v: float(f"{v:.3g}") if v == v else np.nan)
    tests["odds_ratio_death_given_removed"] = tests[
        "odds_ratio_death_given_removed"
    ].map(lambda v: float(f"{v:.3g}") if v == v and np.isfinite(v) else np.nan)
    tests.to_csv(OUT / "restriction_selection_tests.csv", index=False)
    print("\n=== Is removal independent of 28-day death? (Fisher exact) ===")
    print(tests.to_string(index=False))


if __name__ == "__main__":
    selection_tests()
