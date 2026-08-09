from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "03_cluster_robustness"
    / "run_mixak_refit.R"
)


def test_refit_uses_dynamic_features_and_correct_random_intercepts():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "random.intercept = rep(TRUE, length(features))" in text
    assert "fixed_effects <- setNames" in text
    assert "random_effects <- setNames" in text


def test_refit_keeps_probabilities_on_native_scale():
    text = SCRIPT.read_text(encoding="utf-8")
    assert 'probability_scale = "native 0-1 (no division)"' in text
    assert 'assignment_config$probability_scaling != "none"' in text
    assert "quant.comp.prob[[\"50%\"]] / 2" not in text


def test_uncertain_label_cannot_collide_with_candidate_cluster_number():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "uncertain_label <- k + 1L" in text
    assert "group_hpd == uncertain_label" in text
