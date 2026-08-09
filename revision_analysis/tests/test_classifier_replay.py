from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "06_classifier_validation"
    / "replay_archived_autogluon.py"
)
SPEC = importlib.util.spec_from_file_location("replay_archived_autogluon", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_multiclass_brier_is_zero_for_perfect_probabilities():
    y = np.array([1, 2, 3])
    prob = np.eye(3)
    assert MODULE.multiclass_brier(y, prob) == 0


def test_aggregate_outputs_do_not_contain_patient_identifiers():
    y = np.array([1, 1, 2, 2, 3, 3])
    prob = np.array(
        [
            [0.8, 0.1, 0.1],
            [0.7, 0.2, 0.1],
            [0.1, 0.8, 0.1],
            [0.2, 0.7, 0.1],
            [0.1, 0.1, 0.8],
            [0.1, 0.2, 0.7],
        ]
    )
    pred = np.argmax(prob, axis=1) + 1
    outputs = MODULE.aggregate_evaluation("model", "cohort", y, pred, prob)
    for frame in outputs:
        assert "stay_id" not in frame.columns


def test_simple_comparator_is_prespecified_and_small():
    assert MODULE.SIMPLE_FEATURES == [
        "creatinine_min",
        "creatinine_max",
        "crea_divide_basecrea_min",
        "crea_divide_basecrea_max",
        "urineoutput_min",
        "urineoutput_mean",
    ]
