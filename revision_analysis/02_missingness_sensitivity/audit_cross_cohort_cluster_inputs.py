#!/usr/bin/env python3
"""Audit the frozen cross-cohort clustering inputs without changing them.

Only aggregate metadata and summary statistics are written. Patient identifiers and
row-level trajectories are intentionally excluded from revision outputs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_ROOT = REPO_ROOT / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
OUTPUT_ROOT = REPO_ROOT / "02_revision_outputs/reports/W2_cross_cohort_inputs"


@dataclass(frozen=True)
class CohortSpec:
    label: str
    matrix: Path
    baseline: Path
    baseline_id: str
    baseline_column: str
    baseline_multiplier_to_mg_dl: float
    renal_features: tuple[str, ...]


COHORTS = {
    "MIMIC-IV": CohortSpec(
        label="MIMIC-IV",
        matrix=DATA_ROOT / "01.MIMICIV_SAKI_trajCluster/df_mixAK_fea4_C3.csv",
        baseline=DATA_ROOT / "00.data_mimic/disease_definition/AKI/df_base_crea.csv",
        baseline_id="stay_id",
        baseline_column="baseline_Scr",
        baseline_multiplier_to_mg_dl=1.0,
        renal_features=("bun", "creatinine", "urineoutput", "crea_divide_basecrea"),
    ),
    "eICU-CRD": CohortSpec(
        label="eICU-CRD",
        matrix=DATA_ROOT / "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv",
        baseline=DATA_ROOT / "00.data_eicu/disease_definition/AKI/df_base_crea.csv",
        baseline_id="stay_id",
        baseline_column="baseline_creatinine",
        baseline_multiplier_to_mg_dl=1.0,
        renal_features=("bun", "creatinine", "urineoutput", "crea_divide_basecrea"),
    ),
    "AUMC": CohortSpec(
        label="AUMC",
        matrix=DATA_ROOT / "02.AUMCdb_SAKI_trajCluster/df_mixAK_fea3_C3_aumc.csv",
        baseline=DATA_ROOT / "00.data_aumc/disease_definition/AKI/baseline_creatinine.csv",
        baseline_id="admissionid",
        baseline_column="baseline_creatinine",
        baseline_multiplier_to_mg_dl=0.01131,
        renal_features=("creatinine", "urineoutput", "crea_divide_basecrea"),
    ),
}

ALL_RENAL_FEATURES = ("bun", "creatinine", "urineoutput", "crea_divide_basecrea")
FEATURE_UNITS = {
    "bun": "mg/dL",
    "creatinine": "mg/dL",
    "urineoutput": "mL per 6-h window",
    "crea_divide_basecrea": "ratio",
}
GROUP_LABELS = {1: "DR", 2: "RR", 3: "PW", "1": "DR", "2": "RR", "3": "PW"}
RATIO_EXACT_FRACTION_THRESHOLD = 0.999
RATIO_MAX_ERROR_THRESHOLD = 0.01


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_ROOT)
    return parser.parse_args()


def _numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def validate_schema(frame: pd.DataFrame, spec: CohortSpec) -> None:
    required = {"stay_id", "time", "groupHPD", *spec.renal_features}
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"{spec.label}: missing required clustering columns: {missing}")
    if frame.duplicated(["stay_id", "time"]).any():
        raise ValueError(f"{spec.label}: duplicate patient-time rows detected")
    labels = set(frame["groupHPD"].dropna().astype(str))
    if not labels.issubset({"1", "2", "3"}):
        raise ValueError(f"{spec.label}: unexpected phenotype labels: {sorted(labels)}")


def baseline_ratio_audit(frame: pd.DataFrame, spec: CohortSpec) -> dict[str, object]:
    baseline = pd.read_csv(spec.baseline)
    required = {spec.baseline_id, spec.baseline_column}
    if not required.issubset(baseline.columns):
        raise ValueError(
            f"{spec.label}: baseline source lacks {sorted(required - set(baseline.columns))}"
        )
    baseline = baseline[[spec.baseline_id, spec.baseline_column]].copy()
    baseline = baseline.rename(columns={spec.baseline_id: "stay_id"})
    baseline["baseline_mg_dl"] = (
        _numeric(baseline[spec.baseline_column]) * spec.baseline_multiplier_to_mg_dl
    )
    baseline = baseline.dropna(subset=["stay_id", "baseline_mg_dl"])
    if baseline.duplicated("stay_id").any():
        raise ValueError(f"{spec.label}: duplicate identifiers in baseline source")

    linked = frame[["stay_id", "creatinine", "crea_divide_basecrea"]].merge(
        baseline[["stay_id", "baseline_mg_dl"]], on="stay_id", how="left", validate="many_to_one"
    )
    expected = (_numeric(linked["creatinine"]) / linked["baseline_mg_dl"]).round(2)
    observed = _numeric(linked["crea_divide_basecrea"])
    evaluable = expected.notna() & observed.notna()
    error = (expected[evaluable] - observed[evaluable]).abs()
    in_range = baseline["baseline_mg_dl"].between(0.5, 1.5, inclusive="left")
    matrix_ids = pd.Index(frame["stay_id"].dropna().unique())
    baseline_ids = pd.Index(baseline.loc[in_range, "stay_id"].unique())
    return {
        "cohort": spec.label,
        "baseline_source_sha256": sha256(spec.baseline),
        "baseline_source_unit": (
            "micromol/L" if spec.baseline_multiplier_to_mg_dl != 1.0 else "mg/dL"
        ),
        "conversion_to_mg_dl": spec.baseline_multiplier_to_mg_dl,
        "baseline_operational_rule": (
            "single archived baseline field; patients with missing baseline or "
            "baseline outside 0.5 to <1.5 mg/dL were excluded"
        ),
        "matrix_patients": int(len(matrix_ids)),
        "matrix_patients_linked_to_in_range_baseline": int(matrix_ids.isin(baseline_ids).sum()),
        "matrix_rows_evaluable_for_ratio_reconstruction": int(evaluable.sum()),
        "ratio_exact_after_rounding_fraction": (
            float(np.isclose(error, 0.0, atol=1e-12).mean()) if len(error) else np.nan
        ),
        "ratio_absolute_error_max": float(error.max()) if len(error) else np.nan,
        "alternate_baseline_definition_available_in_cluster_archive": False,
    }


def audit_cohort(spec: CohortSpec) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    frame = pd.read_csv(spec.matrix)
    validate_schema(frame, spec)
    time = _numeric(frame["time"])
    patient_counts = frame.groupby("stay_id").size()
    labels = frame[["stay_id", "groupHPD"]].drop_duplicates()
    if labels.duplicated("stay_id").any():
        raise ValueError(f"{spec.label}: patient has more than one phenotype label")

    inventory = {
        "cohort": spec.label,
        "matrix_sha256": sha256(spec.matrix),
        "rows": int(len(frame)),
        "patients": int(frame["stay_id"].nunique()),
        "time_min": float(time.min()),
        "time_max": float(time.max()),
        "unique_time_values": int(time.nunique()),
        "patient_time_duplicates": int(frame.duplicated(["stay_id", "time"]).sum()),
        "minimum_rows_per_patient": int(patient_counts.min()),
        "median_rows_per_patient": float(patient_counts.median()),
        "maximum_rows_per_patient": int(patient_counts.max()),
        "renal_feature_count": len(spec.renal_features),
        "renal_features": ";".join(spec.renal_features),
        "BUN_in_executable_matrix": "bun" in spec.renal_features,
    }

    feature_rows: list[dict[str, object]] = []
    for feature in ALL_RENAL_FEATURES:
        available = feature in frame.columns and feature in spec.renal_features
        values = _numeric(frame[feature]) if feature in frame.columns else pd.Series(dtype=float)
        observed = values.dropna()
        feature_rows.append(
            {
                "cohort": spec.label,
                "feature": feature,
                "unit": FEATURE_UNITS[feature],
                "used_in_executable_clustering_input": available,
                "missing_fraction": float(values.isna().mean()) if len(values) else 1.0,
                "minimum": float(observed.min()) if len(observed) else np.nan,
                "p01": float(observed.quantile(0.01)) if len(observed) else np.nan,
                "median": float(observed.median()) if len(observed) else np.nan,
                "p99": float(observed.quantile(0.99)) if len(observed) else np.nan,
                "maximum": float(observed.max()) if len(observed) else np.nan,
            }
        )

    label_rows = []
    for raw_label, count in labels["groupHPD"].astype(str).value_counts().sort_index().items():
        label_rows.append(
            {
                "cohort": spec.label,
                "phenotype": GROUP_LABELS[raw_label],
                "patients": int(count),
                "fraction": float(count / len(labels)),
            }
        )
    return inventory, feature_rows, label_rows, baseline_ratio_audit(frame, spec)


def determine_audit_status(
    inventory: pd.DataFrame, baselines: pd.DataFrame
) -> dict[str, object]:
    """Derive the release decision from prespecified integrity thresholds."""
    failures: list[str] = []
    for _, row in inventory.iterrows():
        if int(row["patient_time_duplicates"]) != 0:
            failures.append(f"{row['cohort']}: duplicate patient-time rows")
    for _, row in baselines.iterrows():
        cohort = row["cohort"]
        if int(row["matrix_patients_linked_to_in_range_baseline"]) != int(
            row["matrix_patients"]
        ):
            failures.append(f"{cohort}: incomplete linkage to eligible baseline SCr")
        exact_fraction = float(row["ratio_exact_after_rounding_fraction"])
        max_error = float(row["ratio_absolute_error_max"])
        if not np.isfinite(exact_fraction) or exact_fraction < RATIO_EXACT_FRACTION_THRESHOLD:
            failures.append(f"{cohort}: creatinine-ratio exact-match fraction below threshold")
        if not np.isfinite(max_error) or max_error > RATIO_MAX_ERROR_THRESHOLD:
            failures.append(f"{cohort}: creatinine-ratio maximum error above threshold")
    return {
        "overall_status": (
            "FAIL_CROSS_COHORT_INPUT_INTEGRITY"
            if failures
            else "PASS_WITH_DISCLOSED_STRUCTURAL_HETEROGENEITY"
        ),
        "passed": not failures,
        "failure_reasons": failures,
        "thresholds": {
            "patient_time_duplicates": 0,
            "baseline_linkage_fraction": 1.0,
            "ratio_exact_after_rounding_fraction_minimum": RATIO_EXACT_FRACTION_THRESHOLD,
            "ratio_absolute_error_maximum": RATIO_MAX_ERROR_THRESHOLD,
        },
        "cross_cohort_difference": "AUMC lacks BUN and uses three renal variables",
        "baseline_alternative_available": False,
        "patient_level_output_written": False,
    }


def write_report(
    output_dir: Path,
    inventory: pd.DataFrame,
    features: pd.DataFrame,
    baselines: pd.DataFrame,
    status: dict[str, object],
) -> None:
    feature_use = features.pivot(
        index="feature", columns="cohort", values="used_in_executable_clustering_input"
    ).reset_index()
    baseline_display = baselines[
        [
            "cohort",
            "baseline_source_unit",
            "conversion_to_mg_dl",
            "matrix_patients",
            "matrix_patients_linked_to_in_range_baseline",
            "ratio_exact_after_rounding_fraction",
            "ratio_absolute_error_max",
        ]
    ]
    report = f"""# W2 跨队列聚类输入、单位与基线肌酐审计

**状态：{status['overall_status']}**

## 审计结论

- 三个冻结聚类矩阵均无重复 patient-time 行，每位患者只有一个归档表型标签。
- MIMIC-IV 与 eICU-CRD 的可执行聚类输入包含 BUN、肌酐、6 小时尿量和肌酐/基线肌酐比值 4 个肾脏变量。
- AUMC 的归档矩阵不含 BUN，实际使用后 3 个变量。因此修订稿必须描述为“两个四变量模型和一个三变量外部验证模型”，不能声称三个数据库使用完全相同的四变量输入。
- 尿量字段的可审计含义为每个 6 小时窗口内的尿量体积（mL）；不是 24 小时总量，也没有按体重归一化。
- 三队列均使用单一归档基线肌酐字段，入模前要求 0.5–<1.5 mg/dL；缺失或超出范围者被排除，而不是插补。AUMC 原始单位为 µmol/L，按 0.01131 转换为 mg/dL。
- 冻结聚类归档中没有第二套可独立重建的基线肌酐定义。因此本次修订可澄清并审计既有定义，但不能伪造“替代基线定义”敏感性；这项限制应在回复和正文中明确承认。
- 肌酐/基线肌酐比值可由归档肌酐与基线值按两位小数重建，重建统计见下表。

## 矩阵结构

{inventory.to_markdown(index=False)}

## 实际进入聚类的肾脏变量

{feature_use.to_markdown(index=False)}

## 基线肌酐与比值重建

{baseline_display.round(6).to_markdown(index=False)}

## 修订口径

1. Methods 按队列分别列出实际变量，而不是笼统写“所有数据库均为四变量”。
2. 将基线肌酐缺失/范围排除写入纳排流程，并在局限性中说明这限制了对既往肾功能异常人群的外推。
3. 将 AUMC 缺少 BUN 作为跨库结构差异，不把它包装成完全相同模型的独立复制。
4. 单位和数值范围汇总见 `cluster_feature_distribution_audit.csv`；全部结果仅为聚合统计，不含患者标识。
"""
    (output_dir / "W2_CROSS_COHORT_CLUSTER_INPUT_AUDIT.md").write_text(report, encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inventories: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    label_rows: list[dict[str, object]] = []
    baseline_rows: list[dict[str, object]] = []
    for spec in COHORTS.values():
        inventory, cohort_features, cohort_labels, baseline = audit_cohort(spec)
        inventories.append(inventory)
        feature_rows.extend(cohort_features)
        label_rows.extend(cohort_labels)
        baseline_rows.append(baseline)

    inventory_df = pd.DataFrame(inventories)
    features_df = pd.DataFrame(feature_rows)
    labels_df = pd.DataFrame(label_rows)
    baselines_df = pd.DataFrame(baseline_rows)
    inventory_df.to_csv(args.output_dir / "cluster_matrix_inventory.csv", index=False)
    features_df.to_csv(args.output_dir / "cluster_feature_distribution_audit.csv", index=False)
    labels_df.to_csv(args.output_dir / "cluster_label_counts.csv", index=False)
    baselines_df.to_csv(args.output_dir / "baseline_creatinine_ratio_audit.csv", index=False)
    status = determine_audit_status(inventory_df, baselines_df)
    (args.output_dir / "audit_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_report(args.output_dir, inventory_df, features_df, baselines_df, status)
    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 0 if status["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
