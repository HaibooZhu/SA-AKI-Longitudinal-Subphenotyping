#!/usr/bin/env python3
"""Build Supplementary Table S20 -- informative selection by the documentation
density restrictions. Written from the W20 aggregate outputs only."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
REPORT = REPO / "02_revision_outputs/reports/W20_documentation_selection"
OUT = REPO / "02_revision_outputs/tables"

COHORT = {"MIMIC-IV": "MIMIC-IV", "eICU-CRD": "eICU-CRD", "AUMC": "AmsterdamUMCdb"}
RESTRICTION = {
    ">=50% urine-output coverage": "≥50% urine-output coverage",
    "complete 30-window follow-up": "Complete 30-window follow-up",
}
PHENO = {"RR": "Rapid Recovery", "DR": "Delayed Recovery", "PW": "Progressive Worsening", "ALL": "All patients"}


def pct(value: float) -> str:
    return "—" if pd.isna(value) else f"{100 * value:.1f}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    retention = pd.read_csv(REPORT / "restriction_retention.csv")
    tests = pd.read_csv(REPORT / "restriction_selection_tests.csv")
    frame = retention.merge(
        tests[
            [
                "cohort",
                "restriction",
                "phenotype",
                "odds_ratio_death_given_removed",
                "fisher_p",
            ]
        ],
        on=["cohort", "restriction", "phenotype"],
    )
    frame["Cohort"] = frame.cohort.map(COHORT)
    frame["Restriction"] = frame.restriction.map(RESTRICTION)
    frame["Trajectory subphenotype"] = frame.phenotype.map(PHENO)
    table = pd.DataFrame(
        {
            "Cohort": frame["Cohort"],
            "Restriction": frame["Restriction"],
            "Trajectory subphenotype": frame["Trajectory subphenotype"],
            "Patients, n": frame.n_total.astype(int),
            "Retained, n": frame.n_retained.astype(int),
            "Retained, %": frame.retention_rate.map(pct),
            "28-day mortality, removed patients, %": frame.mortality28_removed.map(pct),
            "28-day mortality, retained patients, %": frame.mortality28_retained.map(pct),
            "OR for death given removal": frame.odds_ratio_death_given_removed.map(
                lambda v: "—" if pd.isna(v) else f"{v:.2f}"
            ),
            "P value": frame.fisher_p.map(
                lambda v: "—" if pd.isna(v) else ("<0.001" if v < 0.001 else f"{v:.3f}")
            ),
        }
    )
    order_r = ["≥50% urine-output coverage", "Complete 30-window follow-up"]
    order_c = ["eICU-CRD", "MIMIC-IV", "AmsterdamUMCdb"]
    order_p = list(PHENO.values())
    table = table.assign(
        _r=table["Restriction"].map({v: i for i, v in enumerate(order_r)}),
        _c=table["Cohort"].map({v: i for i, v in enumerate(order_c)}),
        _p=table["Trajectory subphenotype"].map({v: i for i, v in enumerate(order_p)}),
    ).sort_values(["_r", "_c", "_p"]).drop(columns=["_r", "_c", "_p"])

    table.to_csv(OUT / "Table_S20_documentation_selection.csv", index=False)
    caption = (
        "**Table S20. Documentation-density restrictions remove patients "
        "non-randomly with respect to death.** Two sensitivity analyses restrict "
        "the analysed population by documentation density: retaining only eICU-CRD "
        "patients with at least 15 of 30 planned urine-output windows documented, "
        "and retaining only patients observed across all 30 planned windows. For "
        "each restriction the table reports how many patients of each trajectory "
        "subphenotype are retained, and the 28-day mortality of the patients the "
        "restriction removes versus those it retains. Odds ratios and P values are "
        "from Fisher exact tests of the association between removal and 28-day "
        "death within each subphenotype. Patients with missing 28-day mortality "
        "were excluded (eICU-CRD, n = 17). Within Progressive Worsening the removed "
        "patients are consistently those who die, so the restricted populations are "
        "survivor-enriched and cannot test whether that subphenotype exists in the "
        "full cohort. See also Figure S15."
    )
    (OUT / "Table_S20_caption.md").write_text(caption + "\n", encoding="utf-8")
    print(table.to_markdown(index=False))
    print()
    print("wrote", (OUT / "Table_S20_documentation_selection.csv").relative_to(REPO))


if __name__ == "__main__":
    main()
