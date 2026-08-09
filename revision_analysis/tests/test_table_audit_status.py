from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pandas as pd


SCRIPT = Path(__file__).parents[1] / "01_data_audit" / "audit_table_s3.py"
SPEC = importlib.util.spec_from_file_location("audit_table_s3", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def comparison(rows: list[tuple[str, str, str, int]], match: bool = False) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "cohort": cohort,
                "feature": feature,
                "source_group": group,
                "day": day,
                "match": match,
            }
            for cohort, feature, group, day in rows
        ],
        columns=["cohort", "feature", "source_group", "day", "match"],
    )


def passing_comparison() -> pd.DataFrame:
    return comparison([("mimic", "Glucose", "C2", 7)], match=True)


def test_source_mismatch_fails() -> None:
    decision = MODULE.decide_audit(passing_comparison(), comparison([("mimic", "pH", "C2", 7)]))
    assert decision.status == "FAIL"
    assert decision.exit_code != 0


def test_exact_expected_rendering_error_passes_with_qualification() -> None:
    embedded = comparison(sorted(MODULE.EXPECTED_EMBEDDED_MISMATCH_KEYS))
    decision = MODULE.decide_audit(embedded, passing_comparison())
    assert decision.status == "PASS_WITH_EXPECTED_RENDERING_ERROR"
    assert decision.exit_code == 0


def test_unexpected_embedded_mismatch_requires_review() -> None:
    embedded = comparison([("eicu", "Glucose", "C2", 7)])
    decision = MODULE.decide_audit(embedded, passing_comparison())
    assert decision.status == "REVIEW_REQUIRED"
    assert decision.exit_code != 0


def test_no_mismatches_passes() -> None:
    decision = MODULE.decide_audit(passing_comparison(), passing_comparison())
    assert decision.status == "PASS"
    assert decision.exit_code == 0
