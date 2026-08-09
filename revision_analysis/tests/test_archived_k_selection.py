from __future__ import annotations

from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "03_cluster_robustness"
    / "extract_archived_eicu_k_diagnostics.R"
)


def test_archived_k_script_limits_candidates_to_traceable_models():
    text = SCRIPT.read_text(encoding="utf-8")
    assert 'required_objects <- paste0("mod", 2:5)' in text
    assert "mod6" not in text
    assert "mod7" not in text
    assert "mod8" not in text


def test_probability_scaling_is_native_and_patient_outputs_are_absent():
    text = SCRIPT.read_text(encoding="utf-8")
    assert 'probability_scale = "native 0-1 (no division)"' in text
    assert "assignments.csv" not in text
    assert "stay_id" not in text
