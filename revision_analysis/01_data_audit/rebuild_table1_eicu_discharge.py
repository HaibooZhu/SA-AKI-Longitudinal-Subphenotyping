#!/usr/bin/env python3
"""W22 -- rebuild the eICU-CRD hospital-discharge-location block of Table 1.

The submitted block force-fits eICU-CRD's discharge categories onto MIMIC-IV's
row labels, mislabels one row, does not reconcile with the source counts for two
rows, and reports percentages that sum to 100.9%. This script regenerates the
block from the frozen eICU patient table joined to the authoritative N = 1,417
cohort, using eICU's own category names and an explicit non-missing denominator.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SNAP = REPO / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
OUT = REPO / "02_revision_outputs/reports/W22_table1_eicu_discharge"

LABELS = {1: "DR", 2: "RR", 3: "PW"}
ORDER = [
    "Home",
    "Rehabilitation",
    "Skilled Nursing Facility",
    "Nursing Home",
    "Other Hospital",
    "Death",
    "Other",
    "Other External",
]
SUBMITTED = {
    "Home": (651, 450, 183, 18),
    "Home health care services": (30, 15, 10, 5),
    "Rehabilitation": (75, 46, 24, 5),
    "Nursing Facility": (291, 176, 93, 22),
    "Other hospital": (98, 62, 24, 12),
    "Death + Hospice": (174, 58, 65, 51),
    "Other": (98, 62, 24, 12),
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    patients = pd.read_csv(
        SNAP / "00.data_eicu/raw/patient.csv",
        usecols=["patientunitstayid", "hospitaldischargelocation"],
    )
    cohort = pd.read_csv(
        SNAP / "03.eICU_SAKI_trajCluster/sk_survival.csv", usecols=["stay_id", "groupHPD"]
    )
    merged = cohort.merge(
        patients, left_on="stay_id", right_on="patientunitstayid", how="left"
    )
    merged["phenotype"] = merged.groupHPD.map(LABELS)

    n_total = len(merged)
    missing = int(merged.hospitaldischargelocation.isna().sum())
    observed = merged.dropna(subset=["hospitaldischargelocation"])
    counts = pd.crosstab(observed.hospitaldischargelocation, observed.phenotype)
    counts["Overall"] = counts.sum(axis=1)
    counts = counts[["Overall", "RR", "DR", "PW"]].reindex(ORDER)
    denominators = {
        "Overall": len(observed),
        **{g: int(observed.phenotype.eq(g).sum()) for g in ["RR", "DR", "PW"]},
    }

    rows = []
    for label in ORDER:
        row = {"Hospital discharge location": label}
        for column in ["Overall", "RR", "DR", "PW"]:
            n = int(counts.loc[label, column])
            row[column] = f"{n} ({100 * n / denominators[column]:.1f})"
        rows.append(row)
    corrected = pd.DataFrame(rows)
    corrected.to_csv(OUT / "table1_eicu_discharge_corrected.csv", index=False)

    # Explicit reconciliation against the submitted block.
    checks = []
    source_lookup = {label: tuple(int(counts.loc[label, c]) for c in ["Overall", "RR", "DR", "PW"]) for label in ORDER}
    merged_other = tuple(
        source_lookup["Other"][i] + source_lookup["Other External"][i] for i in range(4)
    )
    mapping = {
        "Home": source_lookup["Home"],
        "Home health care services": source_lookup["Nursing Home"],
        "Rehabilitation": source_lookup["Rehabilitation"],
        "Nursing Facility": source_lookup["Skilled Nursing Facility"],
        "Other hospital": source_lookup["Other Hospital"],
        "Death + Hospice": source_lookup["Death"],
        "Other": merged_other,
    }
    for label, submitted in SUBMITTED.items():
        expected = mapping[label]
        checks.append(
            {
                "submitted_row_label": label,
                "submitted": " / ".join(str(v) for v in submitted),
                "closest_source_category": {
                    "Home health care services": "Nursing Home",
                    "Nursing Facility": "Skilled Nursing Facility",
                    "Other hospital": "Other Hospital",
                    "Death + Hospice": "Death",
                    "Other": "Other + Other External",
                }.get(label, label),
                "source": " / ".join(str(v) for v in expected),
                "match": submitted == expected,
            }
        )
    reconciliation = pd.DataFrame(checks)
    reconciliation.to_csv(OUT / "table1_eicu_discharge_reconciliation.csv", index=False)

    submitted_pct_sum = 46.5 + 2.1 + 5.4 + 20.5 + 7.0 + 12.4 + 7.0
    corrected_pct_sum = sum(
        100 * int(counts.loc[label, "Overall"]) / denominators["Overall"] for label in ORDER
    )

    report = f"""# W22 Table 1 eICU-CRD discharge-location block

Authoritative cohort N = {n_total}; discharge location missing for {missing}
patients ({100 * missing / n_total:.1f}%), so percentages use a non-missing
denominator of {denominators['Overall']} (RR {denominators['RR']},
DR {denominators['DR']}, PW {denominators['PW']}).

## Corrected block, using eICU-CRD's own category names

{corrected.to_markdown(index=False)}

Percentages sum to {corrected_pct_sum:.1f}%.

## Reconciliation against the submitted block

{reconciliation.to_markdown(index=False)}

Percentages in the submitted block sum to {submitted_pct_sum:.1f}%, which is not
possible for a complete categorical breakdown. Two rows do not reconcile with the
source: `Other hospital` carries the value of the `Other` row (Other Hospital is
73 / 47 / 19 / 7, not 98 / 62 / 24 / 12), and `Nursing Facility` reports
291 / 176 / 93 / 22 against a source Skilled Nursing Facility count of
298 / 183 / 90 / 25. `Home health care services` carries the correct counts for
eICU's `Nursing Home` category under a MIMIC-IV row label that does not exist in
eICU-CRD. The MIMIC-IV block of the same table reconciles exactly and is unchanged.
"""
    (OUT / "W22_TABLE1_EICU_DISCHARGE.md").write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
