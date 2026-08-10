#!/usr/bin/env python3
"""Regenerate harmonized longitudinal descriptive tables from the audited source matrix."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DAY_MAP = {
    -2: -1, -1: -1,
    1: 1, 2: 1, 3: 1, 4: 1,
    5: 2, 6: 2, 7: 2, 8: 2,
    9: 3, 10: 3, 11: 3, 12: 3,
    13: 4, 14: 4, 15: 4, 16: 4,
    17: 5, 18: 5, 19: 5, 20: 5,
    21: 6, 22: 6, 23: 6, 24: 6,
    25: 7, 26: 7, 27: 7, 28: 7,
}

VARIABLES = [
    ("spo2", "SpO2", "%", "Mean across available 6-h rows"),
    ("fio2", "FiO2", "%", "eICU fraction multiplied by 100; values outside 21-100 removed"),
    ("po2", "PaO2", "mmHg", "Mean across available 6-h rows"),
    ("resp_rate", "Respiratory rate", "breaths/min", "Mean across available 6-h rows"),
    ("sodium", "Sodium", "mmol/L", "Mean across available 6-h rows"),
    ("chloride", "Chloride", "mmol/L", "Mean across available 6-h rows"),
    ("potassium", "Potassium", "mmol/L", "Mean across available 6-h rows"),
    ("calcium", "Calcium", "mg/dL", "Database-derived value retained; cohort-level plausibility audited"),
    ("urineoutput", "Mean 6-h-window urine output", "mL", "Arithmetic mean of 6-h window totals, not a daily sum"),
    ("creatinine", "Creatinine", "mg/dL", "Mean across available 6-h rows"),
    ("crea_divide_basecrea", "Creatinine/baseline creatinine", "ratio", "Mean of the window-level ratio"),
    ("bun", "Blood urea nitrogen", "mg/dL", "Unavailable in AUMC clustering matrix"),
    ("heart_rate", "Heart rate", "beats/min", "Mean across available 6-h rows"),
    ("dbp", "Diastolic blood pressure", "mmHg", "Mean across available 6-h rows"),
    ("sbp", "Systolic blood pressure", "mmHg", "Mean across available 6-h rows"),
    ("mbp", "Mean blood pressure", "mmHg", "Mean across available 6-h rows"),
    ("wbc", "White blood cell count", "10^9/L", "Mean across available 6-h rows"),
    ("temperature", "Temperature", "degrees C", "Mean across available 6-h rows"),
    ("pco2", "PaCO2", "mmHg", "Mean across available 6-h rows"),
    ("baseexcess", "Base excess", "mmol/L", "Mean across available 6-h rows"),
    ("ph", "pH", "unitless", "Mean across available 6-h rows"),
    ("aniongap", "Anion gap", "mmol/L", "Mean across available 6-h rows"),
    ("lactate", "Lactate", "mmol/L", "Mean across available 6-h rows"),
    ("bicarbonate", "Bicarbonate", "mmol/L", "Mean across available 6-h rows"),
    ("glucose", "Glucose", "mg/dL", "Mean across available 6-h rows"),
    ("hematocrit", "Hematocrit", "%", "Fractions 0.05-0.80 multiplied by 100; percentages 5-80 retained; other values removed"),
    ("hemoglobin", "Hemoglobin", "g/dL", "Mean across available 6-h rows"),
    ("platelets", "Platelet count", "10^9/L", "Unavailable in AUMC merged table"),
    ("bilirubin", "Total bilirubin", "mg/dL", "AUMC micromol/L divided by 17.1"),
    ("inr", "International normalized ratio", "ratio", "Unavailable in AUMC/eICU where absent"),
    ("pt", "Prothrombin time", "s", "Mean across available 6-h rows"),
    ("ptt", "Partial thromboplastin time", "s", "Mean across available 6-h rows"),
    ("alp", "Alkaline phosphatase", "U/L", "Unavailable in AUMC merged table"),
    ("ast", "Aspartate aminotransferase", "U/L", "Unavailable in AUMC merged table"),
    ("alt", "Alanine aminotransferase", "U/L", "Unavailable in AUMC merged table"),
]

FEATURES = [item[0] for item in VARIABLES]
DISPLAY = {item[0]: item[1] for item in VARIABLES}
UNITS = {item[0]: item[2] for item in VARIABLES}
TRANSFORM = {item[0]: item[3] for item in VARIABLES}
GROUP = {1: "DR", 2: "RR", 3: "PW"}
COHORT_LABEL = {"mimic": "MIMIC-IV", "aumcdb": "AUMC", "eicu": "eICU-CRD"}


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-matrix",
        type=Path,
        default=(
            repo
            / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
            / "04.other_feature_in_three_dataset/00.data_merge/df_saki_timeseries_feature_all.csv"
        ),
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=repo / "02_revision_outputs/tables/harmonized_longitudinal_tables",
    )
    return parser.parse_args()


def normalize_hematocrit(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    normalized = pd.Series(np.nan, index=values.index, dtype=float)
    fraction = values.between(0.05, 0.80, inclusive="both")
    percent = values.between(5, 80, inclusive="both")
    normalized.loc[fraction] = values.loc[fraction] * 100
    normalized.loc[percent] = values.loc[percent]
    return normalized


def harmonize(source: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    selected = source.loc[source.time.isin(DAY_MAP)].copy()
    selected["day"] = selected.time.map(DAY_MAP)
    before = selected[["dataset"] + FEATURES].isna().groupby(selected.dataset).mean()

    eicu = selected.dataset == "eicu"
    eicu_fio2 = pd.to_numeric(selected.loc[eicu, "fio2"], errors="coerce")
    if eicu_fio2.dropna().gt(1.0).any():
        raise ValueError("eICU FiO2 contains values above 1.0 before fraction-to-percent conversion")
    selected.loc[eicu, "fio2"] = eicu_fio2 * 100
    selected["fio2"] = pd.to_numeric(selected.fio2, errors="coerce").where(
        pd.to_numeric(selected.fio2, errors="coerce").between(21, 100, inclusive="both")
    )
    aumc = selected.dataset == "aumcdb"
    selected.loc[aumc, "bilirubin"] = (
        pd.to_numeric(selected.loc[aumc, "bilirubin"], errors="coerce") / 17.1
    )
    selected["hematocrit"] = normalize_hematocrit(selected.hematocrit)

    after = selected[["dataset"] + FEATURES].isna().groupby(selected.dataset).mean()
    audit = (
        before.stack().rename("missing_fraction_before")
        .to_frame()
        .join(after.stack().rename("missing_fraction_after"))
        .reset_index(names=["cohort", "feature"])
    )
    audit["additional_removed_fraction"] = audit.missing_fraction_after - audit.missing_fraction_before
    return selected, audit


def build_tables(source: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    grouped = source.groupby(["dataset", "groupHPD", "day"], as_index=False)[FEATURES].mean()
    long = grouped.melt(
        id_vars=["dataset", "groupHPD", "day"],
        value_vars=FEATURES,
        var_name="feature",
        value_name="mean",
    )
    long["cohort"] = long.dataset.map(COHORT_LABEL)
    long["phenotype"] = long.groupHPD.map(GROUP)
    long["display_name"] = long.feature.map(DISPLAY)
    long["unit"] = long.feature.map(UNITS)
    long["mean"] = long["mean"].round(4)
    long = long[["cohort", "phenotype", "day", "feature", "display_name", "unit", "mean"]]

    tables = {}
    for cohort in ["MIMIC-IV", "AUMC", "eICU-CRD"]:
        subset = long.loc[long.cohort == cohort].copy()
        wide = subset.pivot(index=["display_name", "unit"], columns=["phenotype", "day"], values="mean")
        desired = [(group, day) for group in ["RR", "DR", "PW"] for day in [-1, 1, 2, 3, 4, 5, 6, 7]]
        wide = wide.reindex(columns=pd.MultiIndex.from_tuples(desired, names=["phenotype", "day"]))
        display_order = [DISPLAY[feature] for feature in FEATURES]
        wide = wide.reindex(display_order, level="display_name")
        tables[cohort] = wide
    return long, tables


def write_report(out_dir: Path, audit: pd.DataFrame) -> None:
    dictionary = pd.DataFrame(
        [
            {
                "source_variable": feature,
                "display_name": DISPLAY[feature],
                "common_unit": UNITS[feature],
                "revision_transform_or_summary": TRANSFORM[feature],
            }
            for feature in FEATURES
        ]
    )
    dictionary.to_csv(out_dir / "longitudinal_variable_dictionary.csv", index=False)
    focus = audit.loc[audit.feature.isin(["fio2", "bilirubin", "hematocrit"])]
    text = f"""# Regenerated longitudinal supplementary tables

## Result

Tables for MIMIC-IV, AUMC, and eICU-CRD were regenerated directly from the audited
source matrix. The MIMIC-IV RR day-7 values now come from source data rather than the
corrupted embedded spreadsheet cells.

## Prespecified harmonization

- eICU FiO2 was converted from fraction to percent and values outside 21-100% were set
  to missing for these descriptive tables.
- AUMC bilirubin was converted from micromol/L to mg/dL by division by 17.1.
- Hematocrit fractions (0.05-0.80) were converted to percent; plausible percentage
  values (5-80) were retained; other values were set to missing. This repairs the mixed
  fraction/percent and outlier contamination found in the AUMC-derived matrix.
- Urine output is explicitly labeled as the arithmetic mean of six-hour window totals.
  It is not labeled as a daily urine-output sum.

## Missingness effect of the unit/range rules

{focus.round(4).to_markdown(index=False)}

## Files

- `harmonized_longitudinal_long.csv`: tidy source data for all tables.
- `Table_S3_MIMIC-IV.csv`, `Table_S4_AUMC.csv`, and `Table_S5_eICU-CRD.csv`: display-ready wide tables.
- `longitudinal_variable_dictionary.csv`: unit and transformation dictionary.
- `harmonization_missingness_audit.csv`: missingness before and after correction rules.
"""
    (out_dir / "REGENERATED_TABLES_README.md").write_text(text, encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(args.source_matrix)
    required_identifiers = {"stay_id", "time", "dataset", "groupHPD"}
    missing_identifiers = sorted(required_identifiers - set(source.columns))
    if missing_identifiers:
        raise ValueError(f"Source matrix is missing identifiers: {missing_identifiers}")
    if source.duplicated(["dataset", "stay_id", "time"]).any():
        raise ValueError("Source matrix contains duplicate cohort-patient-time rows")
    missing = [feature for feature in FEATURES if feature not in source.columns]
    if missing:
        raise ValueError(f"Source matrix is missing expected variables: {missing}")
    harmonized, audit = harmonize(source)
    long, tables = build_tables(harmonized)
    long.to_csv(args.out_dir / "harmonized_longitudinal_long.csv", index=False)
    audit.to_csv(args.out_dir / "harmonization_missingness_audit.csv", index=False)
    for cohort, table in tables.items():
        table.to_csv(args.out_dir / f"Table_S{ {'MIMIC-IV': 3, 'AUMC': 4, 'eICU-CRD': 5}[cohort] }_{cohort}.csv")
    write_report(args.out_dir, audit)
    print(args.out_dir)


if __name__ == "__main__":
    main()
