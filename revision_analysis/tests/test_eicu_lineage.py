from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "00_data_lineage"
    / "build_eicu_lineage_manifest.py"
)
SPEC = importlib.util.spec_from_file_location("build_eicu_lineage_manifest", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def frame(ids, groups=None, outcomes=None):
    data = {"stay_id": ids, "_cid": [MODULE.canonical_id(value) for value in ids]}
    if groups is not None:
        data["groupHPD"] = groups
    if outcomes is not None:
        data["mortality_28d"] = outcomes
    return pd.DataFrame(data)


def test_canonical_id_normalizes_integer_like_values():
    assert MODULE.canonical_id("123.0") == "123"
    assert MODULE.canonical_id(123) == "123"


def test_relation_classification():
    auth = {"1", "2", "3"}
    assert MODULE.classify_relation(auth, {"1", "2", "3"}) == "MATCH"
    assert MODULE.classify_relation(auth, {"1", "2", "3", "4"}) == "SUPERSET_OF_AUTHORITATIVE"
    assert MODULE.classify_relation(auth, {"1", "2"}) == "SUBSET_OF_AUTHORITATIVE"
    assert MODULE.classify_relation(auth, {"2", "4"}) == "PARTIAL_OVERLAP"


def test_comparison_counts_label_and_outcome_mismatch():
    auth = frame([1, 2, 3], [1, 2, 3], [0, 0, 1])
    artifact = frame([1, 2, 3, 4], [1, 3, 3, 1], [0, 1, 1, 0])
    result = MODULE.compare_to_authoritative(auth, artifact)
    assert result["relation"] == "SUPERSET_OF_AUTHORITATIVE"
    assert result["artifact_only_n"] == 1
    assert result["group_mismatch_n"] == 1
    assert result["mortality_mismatch_n"] == 1


def test_final_input_requires_exact_set_and_label_match():
    comparison = {
        "relation": "MATCH",
        "group_mismatch_n": 0,
        "group_conflict_patients": 0,
        "mortality_mismatch_n": 0,
        "mortality_conflict_patients": 0,
    }
    assert MODULE.trust_decision("final_analysis_input", comparison) == "PASS"
    comparison["group_mismatch_n"] = 1
    assert MODULE.trust_decision("final_analysis_input", comparison) == "FAIL"


def test_classifier_frame_reports_ambiguous_and_unknown_ids(tmp_path, monkeypatch):
    model_dir = (
        tmp_path
        / "07.autogluon/01.model"
        / "Result-a1234_selfv2_MimiceICU_AUMC_CorrMICfilt/input"
    )
    model_dir.mkdir(parents=True)
    pd.DataFrame({"stay_id": [1, 3]}).to_csv(model_dir / "train_set.csv", index=False)
    pd.DataFrame({"stay_id": [2]}).to_csv(model_dir / "test_set1.csv", index=False)
    features = pd.DataFrame(
        {
            "stay_id": [1, 1, 2],
            "_cid": ["1", "1", "2"],
            "dataset": ["eicu", "mimic", "eicu"],
        }
    )
    eicu, _, diagnostics = MODULE.build_classifier_eicu_frame(tmp_path, features)
    assert diagnostics["ambiguous_feature_ids"] == 1
    assert diagnostics["unknown_internal_split_ids"] == 1
    assert eicu["_cid"].nunique() == 2
