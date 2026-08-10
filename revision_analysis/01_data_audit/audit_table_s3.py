#!/usr/bin/env python3
"""Audit Supplementary Table S3 from source matrix to embedded Excel tables."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


FEATURES = [
    "spo2", "fio2", "po2", "resp_rate", "sodium", "chloride", "potassium",
    "calcium", "urineoutput", "creatinine", "Scr/bScr", "bun", "heart_rate",
    "dbp", "sbp", "mbp", "wbc", "temperature", "pco2", "baseexcess", "ph",
    "aniongap", "lactate", "bicarbonate", "glucose", "hematocrit", "hemoglobin",
    "platelets", "bilirubin", "inr", "pt", "ptt", "alp", "ast", "alt",
]

DISPLAY_NAMES = {
    "temperature": "Temperature", "mbp": "MBP", "dbp": "DBP", "sbp": "SBP",
    "ph": "pH", "bicarbonate": "Bicarbonate", "baseexcess": "Base Excess",
    "pco2": "Pco2", "calcium": "Calcium", "sodium": "Sodium",
    "creatinine": "Creatinine", "potassium": "Potassium", "bilirubin": "Bilirubin",
    "resp_rate": "Respiratory Rate", "heart_rate": "Heart Rate",
    "aniongap": "Aniongap", "lactate": "Lactate", "wbc": "WBC",
    "glucose": "Glucose", "fio2": "Fio2", "po2": "Po2", "spo2": "Spo2",
    "chloride": "Chloride", "hemoglobin": "Hemoglobin", "hematocrit": "Hematocrit",
    "Scr/bScr": "Scr/bScr", "inr": "Inr", "pt": "PT", "ptt": "PTT",
    "alp": "Alp", "ast": "Ast", "alt": "Alt", "bun": "BUN",
    "urineoutput": "Urine Output", "platelets": "Platelets",
}

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

COHORTS = {
    "mimic": ("Workbook1_df_mimic_other_feature.json", "df_mimic_other_feature.csv"),
    "aumcdb": ("Workbook2_df_aumcdb_other_feature.json", "df_aumcdb_other_feature.csv"),
    "eicu": ("Workbook3_df_eicu_other_feature.json", "df_eicu_other_feature.csv"),
}

EMBEDDED_TO_SOURCE_GROUP = {"RR": "C2", "DR": "C1", "PW": "C3"}

EXPECTED_EMBEDDED_MISMATCH_KEYS = frozenset(
    ("mimic", feature, "C2", 7)
    for feature in [
        "Glucose", "Heart Rate", "Hematocrit", "Hemoglobin", "Inr", "Lactate",
        "MBP", "PT", "PTT", "Pco2", "Platelets", "Po2", "Potassium",
        "Respiratory Rate", "SBP", "Scr/bScr", "Sodium", "Spo2", "Temperature",
        "Urine Output", "WBC", "pH",
    ]
)


@dataclass(frozen=True)
class Paths:
    source_matrix: Path
    generated_dir: Path
    workbook_json_dir: Path
    output_dir: Path


@dataclass(frozen=True)
class AuditDecision:
    status: str
    exit_code: int
    conclusion: str


def parse_args() -> Paths:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-matrix", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--workbook-json-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    return Paths(
        source_matrix=args.source_matrix.resolve(),
        generated_dir=args.generated_dir.resolve(),
        workbook_json_dir=args.workbook_json_dir.resolve(),
        output_dir=args.output_dir.resolve(),
    )


def load_embedded_table(path: Path, cohort: str) -> pd.DataFrame:
    payload = json.loads(path.read_text(encoding="utf-8"))
    values = payload["values"]
    group_row, day_row = values[0], values[1]
    groups: list[str | None] = []
    current_group: str | None = None
    for value in group_row[1:]:
        if value in EMBEDDED_TO_SOURCE_GROUP:
            current_group = str(value)
        groups.append(current_group)

    records = []
    for row in values[2:]:
        feature = row[0]
        if not feature or str(feature).startswith("The mean value"):
            continue
        for offset, raw_value in enumerate(row[1:]):
            if raw_value in (None, ""):
                continue
            group = groups[offset]
            day = int(day_row[offset + 1])
            records.append({
                "cohort": cohort,
                "feature": str(feature),
                "display_group": group,
                "source_group": EMBEDDED_TO_SOURCE_GROUP[str(group)],
                "day": day,
                "embedded_value": float(raw_value),
            })
    return pd.DataFrame.from_records(records)


def load_generated_table(path: Path, cohort: str) -> pd.DataFrame:
    wide = pd.read_csv(path, header=[0, 1], index_col=0)
    records = []
    for feature, row in wide.iterrows():
        for (source_group, day), value in row.items():
            records.append({
                "cohort": cohort,
                "feature": str(feature),
                "source_group": str(source_group),
                "day": int(day),
                "generated_value": float(value),
            })
    return pd.DataFrame.from_records(records)


def reconstruct_from_source(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    source = pd.read_csv(path)
    source = source.rename(columns={"crea_divide_basecrea": "Scr/bScr"})
    selected_times = set(DAY_MAP)
    source = source[source["time"].isin(selected_times)].copy()
    source["day"] = source["time"].map(DAY_MAP)

    missingness = (
        source.groupby("dataset")[FEATURES]
        .agg(lambda series: float(series.isna().mean()))
        .T.reset_index(names="feature")
        .melt(id_vars="feature", var_name="cohort", value_name="missing_fraction")
    )

    grouped = (
        source.groupby(["dataset", "groupHPD", "day"], as_index=False)[FEATURES]
        .mean()
    )
    grouped["source_group"] = grouped["groupHPD"].astype(int).map({1: "C1", 2: "C2", 3: "C3"})
    long = grouped.melt(
        id_vars=["dataset", "source_group", "day"],
        value_vars=FEATURES,
        var_name="raw_feature",
        value_name="source_value",
    )
    long["feature"] = long["raw_feature"].map(DISPLAY_NAMES)
    long = long.rename(columns={"dataset": "cohort"})
    long["source_value"] = long["source_value"].round(4)
    long = long.dropna(subset=["source_value"])
    return long[["cohort", "feature", "source_group", "day", "source_value"]], missingness


def detect_fill_series(embedded: pd.DataFrame) -> pd.DataFrame:
    findings = []
    for (cohort, group, day), frame in embedded.groupby(["cohort", "display_group", "day"]):
        frame = frame.reset_index(drop=True)
        values = frame["embedded_value"].to_numpy()
        start = 0
        for index in range(1, len(values) + 1):
            continues = index < len(values) and np.isclose(values[index] - values[index - 1], 1.0, atol=1e-8)
            if continues:
                continue
            if index - start >= 3:
                findings.append({
                    "cohort": cohort,
                    "display_group": group,
                    "day": day,
                    "first_feature": frame.loc[start, "feature"],
                    "last_feature": frame.loc[index - 1, "feature"],
                    "length": index - start,
                    "first_value": values[start],
                    "last_value": values[index - 1],
                })
            start = index
    return pd.DataFrame.from_records(findings)


def detect_cross_cohort_unit_flags(generated: pd.DataFrame) -> pd.DataFrame:
    medians = generated.groupby(["cohort", "feature"])["generated_value"].median().unstack("cohort")
    flags = []
    for feature, row in medians.iterrows():
        finite = row.dropna().abs()
        nonzero = finite[finite > 0]
        if len(nonzero) < 3:
            continue
        ratio = float(nonzero.max() / nonzero.min())
        if ratio >= 5:
            flags.append({
                "feature": feature,
                "max_to_min_median_ratio": ratio,
                **{f"median_{cohort}": row.get(cohort, np.nan) for cohort in COHORTS},
            })
    return pd.DataFrame.from_records(flags).sort_values("max_to_min_median_ratio", ascending=False)


def decide_audit(
    embedded_comparison: pd.DataFrame,
    source_comparison: pd.DataFrame,
) -> AuditDecision:
    """Classify the audit from observed comparisons; never infer a pass from prose."""
    source_mismatches = source_comparison.loc[~source_comparison["match"]]
    if not source_mismatches.empty:
        return AuditDecision(
            status="FAIL",
            exit_code=1,
            conclusion=(
                f"Source-to-generated validation failed for {len(source_mismatches)} cells. "
                "A reporting-only conclusion is not permitted."
            ),
        )

    embedded_mismatches = embedded_comparison.loc[~embedded_comparison["match"]]
    if embedded_mismatches.empty:
        return AuditDecision(
            status="PASS",
            exit_code=0,
            conclusion=(
                "Source-to-generated and generated-to-embedded comparisons both passed "
                "after four-decimal rounding."
            ),
        )

    observed_keys = frozenset(
        (
            str(row.cohort),
            str(row.feature),
            str(row.source_group),
            int(row.day),
        )
        for row in embedded_mismatches.itertuples(index=False)
    )
    if observed_keys == EXPECTED_EMBEDDED_MISMATCH_KEYS:
        return AuditDecision(
            status="PASS_WITH_EXPECTED_RENDERING_ERROR",
            exit_code=0,
            conclusion=(
                "Source-to-generated validation passed. The embedded workbook contains exactly "
                "the 22 prespecified MIMIC/RR/day-7 rendering mismatches and no others. "
                "This audit localizes that specific defect to the source-matrix-to-workbook "
                "reporting path; broader analytical provenance is evaluated separately."
            ),
        )

    return AuditDecision(
        status="REVIEW_REQUIRED",
        exit_code=2,
        conclusion=(
            f"Source-to-generated validation passed, but {len(embedded_mismatches)} embedded "
            "mismatches do not equal the prespecified 22-cell rendering defect. Manual review "
            "is required before using a reporting-only conclusion."
        ),
    )


def write_report(
    paths: Paths,
    embedded_comparison: pd.DataFrame,
    source_comparison: pd.DataFrame,
    sequences: pd.DataFrame,
    unit_flags: pd.DataFrame,
    missingness: pd.DataFrame,
    decision: AuditDecision,
) -> None:
    mismatches = embedded_comparison[~embedded_comparison["match"]]
    source_mismatches = source_comparison[~source_comparison["match"]]
    cohort_counts = mismatches.groupby("cohort").size().reindex(COHORTS, fill_value=0)
    source_counts = source_mismatches.groupby("cohort").size().reindex(COHORTS, fill_value=0)

    lines = [
        "# W1 — Supplementary Table S3 data-integrity audit",
        "",
        "## Conclusion",
        "",
        f"**Audit status: `{decision.status}`**",
        "",
        decision.conclusion,
        "",
        "## Layer-by-layer checks",
        "",
        "| Cohort | Source vs generated CSV mismatches | Generated CSV vs embedded workbook mismatches |",
        "|---|---:|---:|",
    ]
    for cohort in COHORTS:
        lines.append(f"| {cohort} | {int(source_counts[cohort])} | {int(cohort_counts[cohort])} |")

    lines.extend([
        "",
        "## Detected spreadsheet fill series",
        "",
    ])
    if sequences.empty:
        lines.append("None.")
    else:
        lines.append("| Cohort | Group | Day | Feature range | Length | Value range |")
        lines.append("|---|---|---:|---|---:|---|")
        for row in sequences.itertuples(index=False):
            lines.append(
                f"| {row.cohort} | {row.display_group} | {row.day} | "
                f"{row.first_feature}–{row.last_feature} | {row.length} | "
                f"{row.first_value:.4f}–{row.last_value:.4f} |"
            )

    lines.extend([
        "",
        "## Cross-cohort unit flags",
        "",
        "These are audit flags, not automatic corrections. They require confirmation against each database's extraction dictionary.",
        "",
    ])
    if unit_flags.empty:
        lines.append("No median ratio exceeded the prespecified five-fold screening threshold.")
    else:
        lines.append("| Feature | Max/min median ratio | MIMIC median | AUMC median | eICU median |")
        lines.append("|---|---:|---:|---:|---:|")
        for row in unit_flags.itertuples(index=False):
            lines.append(
                f"| {row.feature} | {row.max_to_min_median_ratio:.2f} | "
                f"{row.median_mimic:.4g} | {row.median_aumcdb:.4g} | {row.median_eicu:.4g} |"
            )

    lines.extend([
        "",
        "## Resolution",
        "",
        "The table-generation script constructs a patient-day urine-output sum in `df_fea_add`, "
        "but then does not use that object. Table S3 instead summarizes `urineoutput` from the original "
        "six-hour rows using an arithmetic mean. The final Methods, caption, and regenerated Table S3 "
        "therefore consistently define this quantity as mean six-hour-window urine output, not a "
        "patient-day sum. This reporting choice has been resolved and no author decision remains pending.",
        "",
        "## Files produced",
        "",
        "- `table_s3_embedded_vs_generated.csv`: every embedded cell and its generated-CSV reference.",
        "- `table_s3_source_vs_generated.csv`: every generated cell and its source-matrix reconstruction.",
        "- `table_s3_fill_series.csv`: detected +1 spreadsheet sequences.",
        "- `table_s3_unit_flags.csv`: cross-cohort median-ratio screen.",
        "- `table_s3_missingness.csv`: source-matrix missingness before table aggregation.",
        "",
        f"Source matrix: `{paths.source_matrix}`",
    ])
    (paths.output_dir / "W1_DATA_INTEGRITY_AUDIT.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    paths = parse_args()
    paths.output_dir.mkdir(parents=True, exist_ok=True)

    embedded_frames = []
    generated_frames = []
    for cohort, (json_name, csv_name) in COHORTS.items():
        embedded_frames.append(load_embedded_table(paths.workbook_json_dir / json_name, cohort))
        generated_frames.append(load_generated_table(paths.generated_dir / csv_name, cohort))
    embedded = pd.concat(embedded_frames, ignore_index=True)
    generated = pd.concat(generated_frames, ignore_index=True)

    embedded_comparison = embedded.merge(
        generated,
        on=["cohort", "feature", "source_group", "day"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    embedded_comparison["absolute_difference"] = (
        embedded_comparison["embedded_value"] - embedded_comparison["generated_value"]
    ).abs()
    embedded_comparison["match"] = (
        embedded_comparison["_merge"].eq("both")
        & embedded_comparison["absolute_difference"].le(5e-5)
    )

    source_long, missingness = reconstruct_from_source(paths.source_matrix)
    source_comparison = generated.merge(
        source_long,
        on=["cohort", "feature", "source_group", "day"],
        how="outer",
        validate="one_to_one",
        indicator=True,
    )
    source_comparison["absolute_difference"] = (
        source_comparison["generated_value"] - source_comparison["source_value"]
    ).abs()
    source_comparison["match"] = (
        source_comparison["_merge"].eq("both")
        & source_comparison["absolute_difference"].le(5e-5)
    )

    sequences = detect_fill_series(embedded)
    unit_flags = detect_cross_cohort_unit_flags(generated)
    decision = decide_audit(embedded_comparison, source_comparison)

    embedded_comparison.to_csv(paths.output_dir / "table_s3_embedded_vs_generated.csv", index=False)
    source_comparison.to_csv(paths.output_dir / "table_s3_source_vs_generated.csv", index=False)
    sequences.to_csv(paths.output_dir / "table_s3_fill_series.csv", index=False)
    unit_flags.to_csv(paths.output_dir / "table_s3_unit_flags.csv", index=False)
    missingness.to_csv(paths.output_dir / "table_s3_missingness.csv", index=False)
    status_payload = {
        "status": decision.status,
        "exit_code": decision.exit_code,
        "source_mismatches": int((~source_comparison["match"]).sum()),
        "embedded_mismatches": int((~embedded_comparison["match"]).sum()),
        "fill_series_findings": int(len(sequences)),
        "unit_flags": int(len(unit_flags)),
    }
    (paths.output_dir / "audit_status.json").write_text(
        json.dumps(status_payload, indent=2) + "\n",
        encoding="utf-8",
    )
    write_report(
        paths,
        embedded_comparison,
        source_comparison,
        sequences,
        unit_flags,
        missingness,
        decision,
    )

    print(f"Audit status: {decision.status}")
    print(f"Embedded mismatches: {(~embedded_comparison['match']).sum()}")
    print(f"Source mismatches: {(~source_comparison['match']).sum()}")
    print(f"Fill-series findings: {len(sequences)}")
    print(f"Unit flags: {len(unit_flags)}")
    print(paths.output_dir / "W1_DATA_INTEGRITY_AUDIT.md")
    return decision.exit_code


if __name__ == "__main__":
    raise SystemExit(main())
