#!/usr/bin/env python3
"""Build revision-only cross-cohort inputs for Editor Concern #2.

Frozen matrices are never modified. Patient-level scenario inputs are written only
to the ignored local output workspace; the tracked/public deliverable is aggregate
code plus an aggregate scenario manifest.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
SNAPSHOT = REPO_ROOT / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
PLANNED_WINDOWS = (-2, -1, *range(1, 29))
COHORTS = {
    "mimic": {
        "label": "MIMIC-IV",
        "matrix": "01.MIMICIV_SAKI_trajCluster/df_mixAK_fea4_C3.csv",
        "ffill": "01.MIMICIV_SAKI_trajCluster/df_im_By_ffill.csv",
        "features": ["bun", "creatinine", "urineoutput", "crea_divide_basecrea"],
        "ffill_columns": {
            "bun_mean": "bun",
            "creatinine_mean": "creatinine",
            "urineoutput_mean": "urineoutput",
        },
        "baseline": "00.data_mimic/disease_definition/AKI/df_base_crea.csv",
        "baseline_id": "stay_id",
        "baseline_column": "baseline_Scr",
        "baseline_multiplier": 1.0,
        "rrt": "00.data_mimic/treatment/lifesupport.csv",
    },
    "eicu": {
        "label": "eICU-CRD",
        "matrix": "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv",
        "ffill": "03.eICU_SAKI_trajCluster/df_im_By_ffill.csv",
        "features": ["bun", "creatinine", "urineoutput", "crea_divide_basecrea"],
        "ffill_columns": {
            "bun": "bun",
            "creatinine": "creatinine",
            "urineoutput": "urineoutput",
        },
        "baseline": "00.data_eicu/disease_definition/AKI/df_base_crea.csv",
        "baseline_id": "stay_id",
        "baseline_column": "baseline_creatinine",
        "baseline_multiplier": 1.0,
        "rrt": "00.data_eicu/treatment/eicu_lifesupport.csv",
    },
    "aumc": {
        "label": "AUMC",
        "matrix": "02.AUMCdb_SAKI_trajCluster/df_mixAK_fea3_C3_aumc.csv",
        "ffill": "02.AUMCdb_SAKI_trajCluster/df_im_By_ffill.csv",
        "features": ["creatinine", "urineoutput", "crea_divide_basecrea"],
        "ffill_columns": {
            "creatinine": "creatinine",
            "urineoutput": "urineoutput",
        },
        "baseline": "00.data_aumc/disease_definition/AKI/baseline_creatinine.csv",
        "baseline_id": "admissionid",
        "baseline_column": "baseline_creatinine",
        "baseline_multiplier": 0.01131,
        "rrt": "00.data_aumc/treatment/aumcdb_lifesupport.csv",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot-root", type=Path, default=SNAPSHOT)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "02_revision_outputs/intermediate/W3_robustness_inputs",
    )
    parser.add_argument(
        "--report-dir",
        type=Path,
        default=REPO_ROOT / "02_revision_outputs/reports/W3_cross_cohort_robustness",
    )
    return parser.parse_args()


def scale_features(
    frame: pd.DataFrame, features: list[str], method: str
) -> tuple[pd.DataFrame, dict[str, float]]:
    scaled = frame.copy()
    scale_values: dict[str, float] = {}
    for feature in features:
        values = pd.to_numeric(scaled[feature], errors="coerce")
        if method == "zscore":
            center = float(values.mean())
            scale = float(values.std(ddof=0))
        elif method == "robust":
            center = float(values.median())
            scale = float(values.quantile(0.75) - values.quantile(0.25))
        else:
            raise ValueError(f"Unknown scaling method: {method}")
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError(f"{feature}: non-positive {method} scale")
        scaled[feature] = (values - center) / scale
        scale_values[f"{feature}_center"] = center
        scale_values[f"{feature}_scale"] = scale
    return scaled, scale_values


def full_followup(frame: pd.DataFrame) -> pd.DataFrame:
    expected = set(PLANNED_WINDOWS)
    eligible = frame.groupby("stay_id")["time"].agg(
        lambda values: len(values) == len(expected) and set(values.astype(int)) == expected
    )
    return frame.loc[frame["stay_id"].isin(eligible[eligible].index)].copy()


def exclude_documented_rrt(
    frame: pd.DataFrame, rrt_path: Path
) -> tuple[pd.DataFrame, int, int]:
    rrt = pd.read_csv(rrt_path, usecols=["stay_id", "is_rrt"])
    if rrt.groupby("stay_id")["is_rrt"].nunique(dropna=True).gt(1).any():
        raise ValueError(f"Conflicting RRT indicators in {rrt_path}")
    rrt = rrt.drop_duplicates("stay_id")
    labels = frame[["stay_id"]].drop_duplicates().merge(
        rrt, on="stay_id", how="left", validate="one_to_one"
    )
    excluded = set(labels.loc[pd.to_numeric(labels.is_rrt, errors="coerce").eq(1), "stay_id"])
    unknown = int(labels.is_rrt.isna().sum())
    return frame.loc[~frame.stay_id.isin(excluded)].copy(), len(excluded), unknown


def limited_forward_fill(
    snapshot: Path, spec: dict[str, object], archived: pd.DataFrame
) -> pd.DataFrame:
    source_columns = ["stay_id", "time", *spec["ffill_columns"].keys()]
    frame = pd.read_csv(snapshot / spec["ffill"], usecols=source_columns).rename(
        columns=spec["ffill_columns"]
    )
    frame = frame.loc[frame.stay_id.isin(archived.stay_id.unique())].copy()
    baseline = pd.read_csv(snapshot / spec["baseline"])[
        [spec["baseline_id"], spec["baseline_column"]]
    ].rename(
        columns={
            spec["baseline_id"]: "stay_id",
            spec["baseline_column"]: "baseline_Scr",
        }
    )
    baseline["baseline_Scr"] = (
        pd.to_numeric(baseline.baseline_Scr, errors="coerce")
        * float(spec["baseline_multiplier"])
    )
    if baseline.duplicated("stay_id").any():
        raise ValueError(f"Duplicate baseline identifiers for {spec['label']}")
    frame = frame.merge(baseline, on="stay_id", how="left", validate="many_to_one")
    frame = frame.loc[frame.baseline_Scr.between(0.5, 1.5, inclusive="left")]
    frame["crea_divide_basecrea"] = (
        pd.to_numeric(frame.creatinine, errors="coerce") / frame.baseline_Scr
    ).round(2)
    labels = archived[["stay_id", "groupHPD"]].drop_duplicates()
    frame = frame.merge(labels, on="stay_id", how="inner", validate="many_to_one")
    frame = frame.loc[pd.to_numeric(frame.time, errors="coerce").isin(PLANNED_WINDOWS)]
    required = ["stay_id", "time", *spec["features"], "groupHPD"]
    frame = frame[required].dropna(subset=spec["features"])
    if frame.duplicated(["stay_id", "time"]).any():
        raise ValueError(f"Duplicate limited-forward-fill rows for {spec['label']}")
    return frame.sort_values(["stay_id", "time"]).reset_index(drop=True)


def manifest_row(
    cohort: str,
    scenario: str,
    frame: pd.DataFrame,
    archived: pd.DataFrame,
    note: str,
) -> dict[str, object]:
    return {
        "cohort": cohort,
        "scenario": scenario,
        "patients": int(frame.stay_id.nunique()),
        "rows": int(len(frame)),
        "retained_patient_fraction": float(
            frame.stay_id.nunique() / archived.stay_id.nunique()
        ),
        "minimum_rows_per_patient": int(frame.groupby("stay_id").size().min()),
        "median_rows_per_patient": float(frame.groupby("stay_id").size().median()),
        "maximum_rows_per_patient": int(frame.groupby("stay_id").size().max()),
        "description": note,
    }


def main() -> int:
    args = parse_args()
    snapshot = args.snapshot_root.resolve()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for cohort, spec in COHORTS.items():
        archived = pd.read_csv(snapshot / spec["matrix"])
        features = list(spec["features"])
        required = ["stay_id", "time", *features, "groupHPD"]
        archived = archived[required].copy()
        scenarios: list[tuple[str, pd.DataFrame, str]] = []
        for method in ["zscore", "robust"]:
            scaled, _ = scale_features(archived, features, method)
            scenarios.append(
                (
                    f"{method}_scaling",
                    scaled,
                    f"Patient-independent global {method} scaling of each renal feature.",
                )
            )
        no_rrt, excluded_n, unknown_n = exclude_documented_rrt(
            archived, snapshot / spec["rrt"]
        )
        scenarios.append(
            (
                "exclude_documented_rrt",
                no_rrt,
                (
                    f"Excluded {excluded_n} patients with documented RRT; retained "
                    f"{unknown_n} patients without an explicit RRT-source row. In this "
                    "exclusion sensitivity, absence of a row is not treated as documented "
                    "RRT and the patient is retained; the unknown count is reported."
                ),
            )
        )
        scenarios.append(
            (
                "complete_30_window_followup",
                full_followup(archived),
                "Restricted to patients containing every one of the 30 planned windows.",
            )
        )
        scenarios.append(
            (
                "limited_forward_fill_complete_rows",
                limited_forward_fill(snapshot, spec, archived),
                (
                    "Used the pre-multiple-imputation forward-fill archive and retained "
                    "only rows complete for all cohort-specific renal features."
                ),
            )
        )
        for scenario, frame, note in scenarios:
            frame.to_csv(args.output_dir / f"{cohort}__{scenario}.csv", index=False)
            rows.append(manifest_row(spec["label"], scenario, frame, archived, note))

    manifest = pd.DataFrame(rows)
    manifest.to_csv(args.report_dir / "cross_cohort_scenario_manifest.csv", index=False)
    report = f"""# Cross-cohort robustness scenario manifest

These inputs were generated in the revision workspace without modifying frozen data.
They operationalize variable scaling, documented RRT exclusion, complete 30-window
follow-up, and a limited-forward-fill/complete-row alternative to multiple imputation.
Patient-level scenario files remain local and are excluded from the public repository.

{manifest.round(3).to_markdown(index=False)}

The RRT scenario excludes documented RRT recipients. Exact first-RRT timestamps are
not available across all three frozen cohorts, so post-RRT truncation cannot be
implemented consistently and is reported as a limitation rather than simulated.
For eICU, absence of an RRT-source row was coded as 0 in the archived outcome notebook;
the present clustering sensitivity uses the narrower estimand "exclude documented
positive RRT" and separately reports absent source rows. This difference in operational
semantics is explicit and is not presented as evidence that absence proves no RRT.
"""
    (args.report_dir / "W3_CROSS_COHORT_SCENARIO_MANIFEST.md").write_text(
        report, encoding="utf-8"
    )
    print(manifest.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
