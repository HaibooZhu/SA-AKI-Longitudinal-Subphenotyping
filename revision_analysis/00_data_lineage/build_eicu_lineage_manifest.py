#!/usr/bin/env python3
"""Build a privacy-preserving eICU artifact lineage manifest.

The script treats the final eICU survival file as the authoritative 1,417-patient
cohort and compares downstream and legacy artifacts against it. Patient IDs are
used only in memory; no patient-level identifiers are written to disk.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Iterable

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SNAPSHOT = (
    REPO_ROOT / "00_frozen_inputs/data_snapshot/remote_project_snapshot"
)
DEFAULT_OUTPUT = REPO_ROOT / "02_revision_outputs/reports/W0_data_lineage"

ID_COLUMN = "stay_id"
COMPARISON_COLUMNS = {
    ID_COLUMN,
    "time",
    "dataset",
    "groupHPD",
    "mortality_28d",
}


@dataclass(frozen=True)
class ArtifactSpec:
    name: str
    role: str
    expected_relation: str
    paths: tuple[Path, ...]
    dataset: str | None = None


def canonical_id(value: object) -> str:
    """Canonicalize integer-like CSV identifiers without exposing them."""
    if pd.isna(value):
        raise ValueError("Missing stay_id encountered")
    text = str(value).strip()
    try:
        number = Decimal(text)
    except InvalidOperation as exc:
        raise ValueError("Non-numeric stay_id encountered") from exc
    if number != number.to_integral_value():
        raise ValueError("Non-integer stay_id encountered")
    return str(int(number))


def file_sha256(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def combined_file_sha256(paths: Iterable[Path]) -> str:
    digest = hashlib.sha256()
    for path in sorted(paths, key=lambda item: item.as_posix()):
        digest.update(path.as_posix().encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_sha256(path).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def id_set_sha256(ids: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for patient_id in sorted(set(ids), key=lambda item: int(item)):
        digest.update(patient_id.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def load_frame(paths: tuple[Path, ...], dataset: str | None = None) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in paths:
        frame = pd.read_csv(path, usecols=lambda column: column in COMPARISON_COLUMNS)
        if ID_COLUMN not in frame.columns:
            raise ValueError(f"Missing {ID_COLUMN} in {path}")
        if dataset is not None:
            if "dataset" not in frame.columns:
                raise ValueError(f"Dataset filter requested but absent in {path}")
            frame = frame.loc[frame["dataset"].eq(dataset)].copy()
        frame["_cid"] = frame[ID_COLUMN].map(canonical_id)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True, sort=False)


def patient_level(frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, int]]:
    conflicts: dict[str, int] = {}
    patient_columns = ["groupHPD", "mortality_28d"]
    for column in patient_columns:
        if column in frame.columns:
            counts = frame.groupby("_cid", sort=False)[column].nunique(dropna=True)
            conflicts[column] = int(counts.gt(1).sum())
    keep = ["_cid"] + [column for column in patient_columns if column in frame.columns]
    patients = frame[keep].drop_duplicates("_cid", keep="first").copy()
    return patients, conflicts


def group_counts(frame: pd.DataFrame) -> str:
    if "groupHPD" not in frame.columns:
        return "{}"
    patients, _ = patient_level(frame)
    values = pd.to_numeric(patients["groupHPD"], errors="coerce")
    counts = values.value_counts(dropna=False).sort_index()
    payload = {
        ("missing" if pd.isna(key) else str(int(key))): int(value)
        for key, value in counts.items()
    }
    return json.dumps(payload, sort_keys=True)


def declared_grain(frame: pd.DataFrame) -> tuple[str, int]:
    if "time" in frame.columns and frame["time"].notna().any():
        grain = ["_cid", "time"]
        return "patient_time", int(frame.duplicated(grain).sum())
    return "patient", int(frame.duplicated(["_cid"]).sum())


def classify_relation(authoritative_ids: set[str], artifact_ids: set[str]) -> str:
    if artifact_ids == authoritative_ids:
        return "MATCH"
    if authoritative_ids < artifact_ids:
        return "SUPERSET_OF_AUTHORITATIVE"
    if artifact_ids < authoritative_ids:
        return "SUBSET_OF_AUTHORITATIVE"
    if authoritative_ids.isdisjoint(artifact_ids):
        return "DISJOINT"
    return "PARTIAL_OVERLAP"


def compare_to_authoritative(
    authoritative: pd.DataFrame, artifact: pd.DataFrame
) -> dict[str, int | str]:
    auth_patients, _ = patient_level(authoritative)
    artifact_patients, conflicts = patient_level(artifact)
    auth_ids = set(auth_patients["_cid"])
    artifact_ids = set(artifact_patients["_cid"])
    overlap = auth_ids & artifact_ids
    result: dict[str, int | str] = {
        "relation": classify_relation(auth_ids, artifact_ids),
        "authoritative_n": len(auth_ids),
        "artifact_n": len(artifact_ids),
        "overlap_n": len(overlap),
        "authoritative_only_n": len(auth_ids - artifact_ids),
        "artifact_only_n": len(artifact_ids - auth_ids),
        "group_conflict_patients": conflicts.get("groupHPD", 0),
        "mortality_conflict_patients": conflicts.get("mortality_28d", 0),
        "group_mismatch_n": 0,
        "mortality_mismatch_n": 0,
    }
    merged = auth_patients.merge(
        artifact_patients,
        on="_cid",
        how="inner",
        suffixes=("_authoritative", "_artifact"),
        validate="one_to_one",
    )
    for column, output_name in [
        ("groupHPD", "group_mismatch_n"),
        ("mortality_28d", "mortality_mismatch_n"),
    ]:
        left = f"{column}_authoritative"
        right = f"{column}_artifact"
        if left not in merged.columns or right not in merged.columns:
            continue
        lhs = pd.to_numeric(merged[left], errors="coerce")
        rhs = pd.to_numeric(merged[right], errors="coerce")
        comparable = lhs.notna() & rhs.notna()
        result[output_name] = int(lhs[comparable].ne(rhs[comparable]).sum())
    return result


def trust_decision(role: str, comparison: dict[str, int | str]) -> str:
    relation = comparison["relation"]
    label_issues = int(comparison["group_mismatch_n"]) + int(
        comparison["group_conflict_patients"]
    )
    outcome_issues = int(comparison["mortality_mismatch_n"]) + int(
        comparison["mortality_conflict_patients"]
    )
    if role == "authoritative_reference":
        return "REFERENCE"
    if role == "final_analysis_input":
        return "PASS" if relation == "MATCH" and label_issues == 0 else "FAIL"
    if role == "final_outcome_input":
        return (
            "PASS"
            if relation == "MATCH" and label_issues == 0 and outcome_issues == 0
            else "FAIL"
        )
    if role == "legacy_or_upstream_snapshot":
        return (
            "NOT_AUTHORIZED_FOR_FINAL_RESULTS"
            if relation != "MATCH" or label_issues or outcome_issues
            else "MATCH_BUT_NOT_AUTHORITATIVE"
        )
    raise ValueError(f"Unknown artifact role: {role}")


def relative_paths(paths: tuple[Path, ...], snapshot: Path) -> str:
    return "; ".join(path.relative_to(snapshot).as_posix() for path in paths)


def make_manifest_row(
    spec: ArtifactSpec,
    frame: pd.DataFrame,
    authoritative: pd.DataFrame,
    snapshot: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    patients, conflicts = patient_level(frame)
    grain, duplicate_rows = declared_grain(frame)
    comparison = compare_to_authoritative(authoritative, frame)
    decision = trust_decision(spec.role, comparison)
    manifest = {
        "artifact": spec.name,
        "role": spec.role,
        "expected_relation": spec.expected_relation,
        "trust_decision": decision,
        "source_paths": relative_paths(spec.paths, snapshot),
        "source_file_count": len(spec.paths),
        "row_count": len(frame),
        "unique_patient_count": len(patients),
        "declared_grain": grain,
        "duplicate_rows_at_declared_grain": duplicate_rows,
        "group_conflict_patients": conflicts.get("groupHPD", 0),
        "mortality_conflict_patients": conflicts.get("mortality_28d", 0),
        "group_counts": group_counts(frame),
        "id_set_sha256": id_set_sha256(patients["_cid"]),
        "source_files_sha256": combined_file_sha256(spec.paths),
    }
    pairwise = {"artifact": spec.name, **comparison, "trust_decision": decision}
    return manifest, pairwise


def build_classifier_eicu_frame(
    snapshot: Path,
    feature_frame: pd.DataFrame,
) -> tuple[pd.DataFrame, tuple[Path, ...], dict[str, int]]:
    model_dir = (
        snapshot
        / "07.autogluon/01.model"
        / "Result-a1234_selfv2_MimiceICU_AUMC_CorrMICfilt/input"
    )
    paths = (model_dir / "train_set.csv", model_dir / "test_set1.csv")
    splits = load_frame(paths)
    feature_membership = feature_frame[["_cid", "dataset"]].drop_duplicates()
    membership_counts = feature_membership.groupby("_cid")["dataset"].nunique()
    ambiguous_ids = set(membership_counts[membership_counts.gt(1)].index)
    eicu_ids = set(feature_membership.loc[feature_membership.dataset.eq("eicu"), "_cid"])
    mimic_ids = set(feature_membership.loc[feature_membership.dataset.eq("mimic"), "_cid"])
    unknown_ids = set(splits["_cid"]) - eicu_ids - mimic_ids
    eicu = splits.loc[splits["_cid"].isin(eicu_ids)].copy()
    diagnostics = {
        "ambiguous_feature_ids": len(ambiguous_ids),
        "unknown_internal_split_ids": len(unknown_ids),
        "eicu_internal_split_patients": eicu["_cid"].nunique(),
    }
    return eicu, paths, diagnostics


def write_markdown(
    output: Path,
    manifest: pd.DataFrame,
    pairwise: pd.DataFrame,
    status: dict[str, object],
) -> None:
    display_columns = [
        "artifact",
        "unique_patient_count",
        "trust_decision",
        "expected_relation",
    ]
    comparison_columns = [
        "artifact",
        "relation",
        "overlap_n",
        "authoritative_only_n",
        "artifact_only_n",
        "group_mismatch_n",
        "mortality_mismatch_n",
    ]
    by_artifact = pairwise.set_index("artifact")

    def comparison_line(artifact: str, label: str) -> str:
        row = by_artifact.loc[artifact]
        return (
            f"- {label}含 {int(row['artifact_n']):,} 名患者；相对权威队列，"
            f"额外 {int(row['artifact_only_n']):,} 人、缺少 "
            f"{int(row['authoritative_only_n']):,} 人，重叠患者中有 "
            f"{int(row['group_mismatch_n']):,} 个表型标签差异和 "
            f"{int(row['mortality_mismatch_n']):,} 个 28 天死亡结局差异。"
        )

    final_inputs = manifest.loc[
        manifest["role"].isin(["final_analysis_input", "final_outcome_input"])
    ]
    passed_final = int(final_inputs["trust_decision"].eq("PASS").sum())
    conclusion_lines = [
        (
            f"- 共核查 {len(final_inputs)} 个最终分析输入，其中 {passed_final} 个通过"
            "患者集合与标签一致性规则。"
        ),
        comparison_line(
            "eicu_legacy_survival_time_summary", "旧版生存时间汇总"
        ),
        comparison_line(
            "eicu_legacy_risk_factor_union", "旧版风险因素合并文件"
        ),
        comparison_line(
            "eicu_legacy_classifier_feature_source", "分类器上游特征快照"
        ),
        (
            "- 分类器内部输入的队列归属检查发现 "
            f"{status['classifier_membership_diagnostics']['ambiguous_feature_ids']} 个"
            "跨数据库歧义 ID 和 "
            f"{status['classifier_membership_diagnostics']['unknown_internal_split_ids']} 个"
            "无法归属的内部拆分 ID；任一计数非零都会使审计失败。"
        ),
        (
            "- 因此，所有标记为 `NOT_AUTHORIZED_FOR_FINAL_RESULTS` 的旧派生文件"
            "只可用于追溯；最终分析必须从权威队列及经核验的上游变量重新派生。"
        ),
    ]

    lines = [
        "# W0 eICU 数据血缘与队列一致性报告",
        "",
        f"**总体状态：{status['overall_status']}**",
        "",
        "本报告将 `03.eICU_SAKI_trajCluster/sk_survival.csv` 中的 1,417 名患者定义为本轮修订的权威 eICU 队列。患者标识只在内存中用于集合比较；所有落盘文件均只含聚合计数和不可逆 SHA-256 摘要。",
        "",
        "## 结论先行",
        "",
        *conclusion_lines,
        "",
        "## 工件清单",
        "",
        manifest[display_columns].to_markdown(index=False),
        "",
        "## 与权威队列的逐项比较",
        "",
        pairwise[comparison_columns].to_markdown(index=False),
        "",
        "## 放行规则",
        "",
        "1. 最终分析输入必须与权威患者集合完全一致。",
        "2. 含表型标签的最终输入不得存在患者内冲突或与权威标签不一致。",
        "3. 含 28 天死亡结局的最终结局输入不得存在患者内冲突或与权威结局不一致。",
        "4. 标记为 `NOT_AUTHORIZED_FOR_FINAL_RESULTS` 的旧快照只可用于追溯问题，不得直接进入修订分析。",
        "",
        "详细文件摘要见 `eicu_artifact_manifest.csv`，逐项集合比较见 `eicu_pairwise_lineage.csv`，机器可读状态见 `eicu_lineage_status.json`。",
        "",
    ]
    (output / "W0_EICU_DATA_LINEAGE.md").write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshot-root", type=Path, default=DEFAULT_SNAPSHOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    snapshot = args.snapshot_root.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    authoritative_spec = ArtifactSpec(
        name="eicu_authoritative_survival_1417",
        role="authoritative_reference",
        expected_relation="REFERENCE",
        paths=(snapshot / "03.eICU_SAKI_trajCluster/sk_survival.csv",),
    )
    authoritative = load_frame(authoritative_spec.paths)
    if authoritative["_cid"].nunique() != 1417:
        raise RuntimeError("Authoritative eICU cohort is not exactly 1,417 patients")

    specs = [
        authoritative_spec,
        ArtifactSpec(
            "eicu_final_cluster_input",
            "final_analysis_input",
            "MATCH",
            (snapshot / "03.eICU_SAKI_trajCluster/df_mixAK_fea4_C3_eicu.csv",),
        ),
        ArtifactSpec(
            "eicu_final_imputed_matrix",
            "final_analysis_input",
            "MATCH",
            (snapshot / "03.eICU_SAKI_trajCluster/df_im_By_MI.csv",),
        ),
        ArtifactSpec(
            "eicu_cross_cohort_longitudinal_source",
            "final_analysis_input",
            "MATCH",
            (
                snapshot
                / "04.other_feature_in_three_dataset/00.data_merge"
                / "df_saki_timeseries_feature_all.csv",
            ),
            dataset="eicu",
        ),
        ArtifactSpec(
            "eicu_legacy_survival_time_summary",
            "legacy_or_upstream_snapshot",
            "TRACE_ONLY",
            (
                snapshot
                / "04.other_feature_in_three_dataset/04.survival_treatment"
                / "df_time_all.csv",
            ),
            dataset="eicu",
        ),
        ArtifactSpec(
            "eicu_legacy_risk_factor_union",
            "legacy_or_upstream_snapshot",
            "TRACE_ONLY",
            tuple(
                snapshot
                / "04.other_feature_in_three_dataset/08.subphenotype_association_analysis"
                / f"df_eicu_c{group}_riskfactor.csv"
                for group in (1, 2, 3)
            ),
        ),
        ArtifactSpec(
            "eicu_legacy_classifier_feature_source",
            "legacy_or_upstream_snapshot",
            "TRACE_ONLY",
            (
                snapshot
                / "07.autogluon/00.data_generate"
                / "df_saki_tsfresh_generate_features_a1234_CorrMICfilt.csv",
            ),
            dataset="eicu",
        ),
    ]

    loaded: dict[str, pd.DataFrame] = {
        authoritative_spec.name: authoritative,
    }
    for spec in specs[1:]:
        loaded[spec.name] = load_frame(spec.paths, dataset=spec.dataset)

    full_feature_path = (
        snapshot
        / "07.autogluon/00.data_generate"
        / "df_saki_tsfresh_generate_features_a1234_CorrMICfilt.csv"
    )
    full_feature = load_frame((full_feature_path,))
    classifier_eicu, classifier_paths, classifier_diagnostics = (
        build_classifier_eicu_frame(snapshot, full_feature)
    )
    classifier_spec = ArtifactSpec(
        "eicu_classifier_actual_internal_inputs",
        "final_analysis_input",
        "MATCH",
        classifier_paths,
    )
    specs.append(classifier_spec)
    loaded[classifier_spec.name] = classifier_eicu

    manifest_rows: list[dict[str, object]] = []
    pairwise_rows: list[dict[str, object]] = []
    for spec in specs:
        manifest, pairwise = make_manifest_row(
            spec, loaded[spec.name], authoritative, snapshot
        )
        manifest_rows.append(manifest)
        pairwise_rows.append(pairwise)

    manifest_df = pd.DataFrame(manifest_rows)
    pairwise_df = pd.DataFrame(pairwise_rows)
    required = manifest_df[manifest_df["role"].isin(["final_analysis_input", "final_outcome_input"])]
    failed_required = required.loc[required["trust_decision"].ne("PASS"), "artifact"].tolist()
    membership_failures = [
        name
        for name, count in classifier_diagnostics.items()
        if name in {"ambiguous_feature_ids", "unknown_internal_split_ids"}
        and int(count) > 0
    ]
    if membership_failures:
        failed_required.append("eicu_classifier_membership_resolution")
    quarantined = manifest_df.loc[
        manifest_df["trust_decision"].eq("NOT_AUTHORIZED_FOR_FINAL_RESULTS"),
        "artifact",
    ].tolist()
    status = {
        "overall_status": "FAIL" if failed_required else (
            "ACTION_REQUIRED" if quarantined else "PASS"
        ),
        "authoritative_patient_count": authoritative["_cid"].nunique(),
        "failed_required_artifacts": failed_required,
        "quarantined_legacy_artifacts": quarantined,
        "classifier_membership_diagnostics": classifier_diagnostics,
        "classifier_membership_failure_reasons": membership_failures,
        "privacy": "No patient-level identifiers are written by this audit.",
        "next_action": (
            "Re-derive outcome and risk-factor analyses from the authoritative cohort; "
            "do not use quarantined 1,970/1,748-patient snapshots."
            if quarantined
            else "All checked artifacts are eligible for their declared role."
        ),
    }

    manifest_df.to_csv(output / "eicu_artifact_manifest.csv", index=False)
    pairwise_df.to_csv(output / "eicu_pairwise_lineage.csv", index=False)
    (output / "eicu_lineage_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    write_markdown(output, manifest_df, pairwise_df, status)

    print(json.dumps(status, indent=2, ensure_ascii=False))
    return 1 if failed_required else 0


if __name__ == "__main__":
    raise SystemExit(main())
