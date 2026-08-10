#!/usr/bin/env python3
"""Publish a whitelisted, identifier-free aggregate revision evidence bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

import pandas as pd


PROHIBITED_COLUMNS = {
    "stay_id",
    "subject_id",
    "hadm_id",
    "admissionid",
    "patientunitstayid",
}
EVIDENCE_FILES = (
    "W0_data_lineage/W0_EICU_DATA_LINEAGE.md",
    "W0_data_lineage/eicu_artifact_manifest.csv",
    "W0_data_lineage/eicu_pairwise_lineage.csv",
    "W0_data_lineage/eicu_lineage_status.json",
    "W1_eicu_uo_sensitivity/eicu_uo_coverage_by_phenotype.csv",
    "W1_eicu_uo_sensitivity/eicu_uo_scenario_status.json",
    "W1_data_integrity/W1_DATA_INTEGRITY_AUDIT.md",
    "W1_data_integrity/audit_status.json",
    "W1_data_integrity/table_s3_embedded_vs_generated.csv",
    "W1_data_integrity/table_s3_fill_series.csv",
    "W1_data_integrity/table_s3_missingness.csv",
    "W1_data_integrity/table_s3_source_vs_generated.csv",
    "W1_data_integrity/table_s3_unit_flags.csv",
    "W2_cross_cohort_inputs/W2_CROSS_COHORT_CLUSTER_INPUT_AUDIT.md",
    "W2_cross_cohort_inputs/cluster_matrix_inventory.csv",
    "W2_cross_cohort_inputs/cluster_feature_distribution_audit.csv",
    "W2_cross_cohort_inputs/cluster_label_counts.csv",
    "W2_cross_cohort_inputs/baseline_creatinine_ratio_audit.csv",
    "W2_cross_cohort_inputs/audit_status.json",
    "W3_archived_k_selection/W3_ARCHIVED_EICU_K_SELECTION.md",
    "W3_archived_k_selection/archived_eicu_k2_k5_diagnostics.csv",
    "W3_archived_k_selection/archived_eicu_k_selection_status.json",
    "W3_fresh_k_grid/W3_FRESH_CROSS_COHORT_K_GRID.md",
    "W3_fresh_k_grid/fresh_k_grid_all_seed_diagnostics.csv",
    "W3_fresh_k_grid/fresh_k_grid_summary.csv",
    "W3_fresh_k_grid/fresh_k_grid_status.json",
    "W3_deep_k_stability/W3_DEEP_K2_K3_STABILITY.md",
    "W3_deep_k_stability/deep_k_fit_diagnostics.csv",
    "W3_deep_k_stability/deep_k_k3_vs_archived.csv",
    "W3_deep_k_stability/deep_k_k3_pairwise_stability.csv",
    "W3_deep_k_stability/deep_k_cohort_summary.csv",
    "W3_deep_k_stability/deep_k_status.json",
    "W3_deep_k_stability/Figure_S2_cross_cohort_k_stability.png",
    "W3_deep_k_stability/Figure_S2_cross_cohort_k_stability.pdf",
    "W3_cross_cohort_robustness/W3_CROSS_COHORT_SCENARIO_MANIFEST.md",
    "W3_cross_cohort_robustness/cross_cohort_scenario_manifest.csv",
    "W3_cross_cohort_robustness/W3_CROSS_COHORT_ROBUSTNESS_REFITS.md",
    "W3_cross_cohort_robustness/cross_cohort_robustness_metrics.csv",
    "W3_cross_cohort_robustness/robustness_refit_status.json",
    "W3_cross_cohort_robustness/cross_cohort_robustness_matrix.png",
    "W3_cross_cohort_robustness/cross_cohort_robustness_matrix.pdf",
    "W3_deep_robustness/W3_DEEP_CAUTION_RERUNS.md",
    "W3_deep_robustness/deep_robustness_seed_metrics.csv",
    "W3_deep_robustness/deep_robustness_scenario_summary.csv",
    "W3_deep_robustness/deep_robustness_status.json",
    "W3_uo_multiseed_sensitivity/W3_UO_MULTI_SEED_SENSITIVITY.md",
    "W3_uo_multiseed_sensitivity/uo_multiseed_seed_metrics.csv",
    "W3_uo_multiseed_sensitivity/uo_multiseed_scenario_summary.csv",
    "W3_uo_multiseed_sensitivity/uo_multiseed_status.json",
    "W3_mixak_provenance/W3_MIXAK_PROVENANCE.md",
    "W3_mixak_provenance/mixak_code_provenance.csv",
    "W3_mixak_provenance/mixak_rdata_provenance.csv",
    "W3_mixak_provenance/mixak_provenance_status.json",
    "W4_classifier_validation/ARCHIVED_AUTOGLUON_REPLAY.md",
    "W4_classifier_validation/primary_and_comparator_metrics.csv",
    "W4_classifier_validation/primary_and_comparator_class_metrics.csv",
    "W4_classifier_validation/primary_and_comparator_pairwise_auc.csv",
    "W4_classifier_validation/primary_and_comparator_calibration.csv",
    "W4_classifier_validation/paired_incremental_value.csv",
    "W4_classifier_validation/archived_replay_environment.json",
    "W5_independent_outcomes/W5_INDEPENDENT_OUTCOMES.md",
    "W5_independent_outcomes/outcome_source_lineage.csv",
    "W5_independent_outcomes/covariate_source_verification.csv",
    "W5_independent_outcomes/outcome_populations.csv",
    "W5_independent_outcomes/outcome_adjusted_effects.csv",
    "W5_independent_outcomes/mortality_standardized_risks.csv",
    "W5_independent_outcomes/outcome_model_missingness.csv",
    "W6_diuretic_exploratory/W6_DIURETIC_EXPLORATORY.md",
    "W6_diuretic_exploratory/historical_matching_spec.csv",
    "W6_diuretic_exploratory/archived_psm_balance_smd.csv",
    "W6_diuretic_exploratory/early_first_dose_descriptive.csv",
    "W6_diuretic_exploratory/early_first_dose_pooled_models.csv",
    "W6_diuretic_exploratory/early_first_dose_interaction_models.csv",
    "W6_diuretic_exploratory/early_first_dose_interaction_tests.csv",
    "W6_diuretic_exploratory/W6_early_diuretic_response_forest.png",
    "W6_diuretic_exploratory/W6_early_diuretic_response_forest.pdf",
    "W8_revision_document_qa/W8_REVISION_DOCUMENT_QA.md",
    "W8_revision_document_qa/revision_document_qa_status.json",
    "W9_numerical_consistency/W9_NUMERICAL_CONSISTENCY_AUDIT.md",
    "W9_numerical_consistency/numerical_consistency_status.json",
    "W10_reviewer_traceability/W10_FINAL_REVIEWER_TRACEABILITY.md",
    "W10_reviewer_traceability/reviewer_traceability_status.json",
)


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--report-root",
        type=Path,
        default=repo / "02_revision_outputs/reports",
    )
    parser.add_argument("--public-output", type=Path, required=True)
    return parser.parse_args()


def audit_file(path: Path) -> dict[str, object]:
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
        prohibited = sorted(PROHIBITED_COLUMNS & set(frame.columns))
        if prohibited:
            raise ValueError(f"{path}: prohibited identifier columns {prohibited}")
        return {"rows": len(frame), "columns": len(frame.columns)}
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        serialized = json.dumps(payload)
        prohibited = sorted(
            column for column in PROHIBITED_COLUMNS if f'"{column}"' in serialized
        )
        if prohibited:
            raise ValueError(f"{path}: prohibited identifier keys {prohibited}")
        return {"json_valid": True}
    if path.suffix.lower() in {".png", ".pdf"}:
        return {
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }
    return {"text_file": True}


def sanitize_text(text: str, report_root: Path) -> str:
    """Remove workstation-specific paths from otherwise public aggregate evidence."""
    revision_root = report_root.resolve().parents[1]
    text = text.replace(str(revision_root), "<REVISION_REPOSITORY>")
    return re.sub(
        r"/(?:Users|Volumes|home)/[^\s`\"']+",
        "<LOCAL_PATH>",
        text,
    )


def publish_file(source: Path, destination: Path, report_root: Path) -> None:
    if source.suffix.lower() in {".md", ".json", ".csv"}:
        text = source.read_text(encoding="utf-8")
        destination.write_text(sanitize_text(text, report_root), encoding="utf-8")
    else:
        shutil.copy2(source, destination)


def main() -> int:
    args = parse_args()
    args.public_output.mkdir(parents=True, exist_ok=True)
    manifest = []
    missing = []
    for relative in EVIDENCE_FILES:
        source = args.report_root / relative
        if not source.exists():
            missing.append(relative)
            continue
        audit = audit_file(source)
        destination = args.public_output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        publish_file(source, destination, args.report_root)
        manifest.append({"path": relative, **audit})
    if missing:
        raise FileNotFoundError(
            "Aggregate evidence bundle is incomplete: " + ", ".join(missing)
        )
    pd.DataFrame(manifest).to_csv(
        args.public_output / "AGGREGATE_EVIDENCE_MANIFEST.csv", index=False
    )
    readme = """# Aggregate revision evidence

This directory contains only whitelisted aggregate outputs generated from the frozen
authorized analysis snapshot. Patient-level inputs, assignments, and predictions are
excluded. `AGGREGATE_EVIDENCE_MANIFEST.csv` records the published file set and basic
shape checks. Re-run the publisher from the exact revision commit to refresh it.
"""
    (args.public_output / "README.md").write_text(readme, encoding="utf-8")
    print(f"Published {len(manifest)} identifier-free aggregate files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
