from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "07_tables_figures"
    / "regenerate_supplementary_tables.py"
)
SPEC = importlib.util.spec_from_file_location("regenerate_supplementary_tables", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def source_fixture() -> pd.DataFrame:
    rows = []
    for dataset, fio2, bilirubin, hematocrit in [
        ("eicu", 0.4, 1.0, 32.0),
        ("aumcdb", 40.0, 17.1, 0.35),
        ("mimic", 50.0, 1.0, 30.0),
    ]:
        row = {feature: 1.0 for feature in MODULE.FEATURES}
        row.update(
            {
                "dataset": dataset,
                "stay_id": len(rows) + 1,
                "time": 1,
                "groupHPD": 1,
                "fio2": fio2,
                "bilirubin": bilirubin,
                "hematocrit": hematocrit,
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def test_cohort_specific_unit_harmonization():
    harmonized, _ = MODULE.harmonize(source_fixture())
    by_cohort = harmonized.set_index("dataset")
    assert by_cohort.loc["eicu", "fio2"] == pytest.approx(40.0)
    assert by_cohort.loc["aumcdb", "bilirubin"] == pytest.approx(1.0)
    assert by_cohort.loc["aumcdb", "hematocrit"] == pytest.approx(35.0)
    assert by_cohort.loc["mimic", "fio2"] == pytest.approx(50.0)


def test_eicu_fio2_above_fraction_scale_is_rejected():
    source = source_fixture()
    source.loc[source.dataset.eq("eicu"), "fio2"] = 40.0
    with pytest.raises(ValueError, match="above 1.0"):
        MODULE.harmonize(source)


def test_hematocrit_outliers_are_removed():
    values = pd.Series([0.35, 35.0, 0.0, 161.0])
    normalized = MODULE.normalize_hematocrit(values)
    assert normalized.iloc[0] == pytest.approx(35.0)
    assert normalized.iloc[1] == pytest.approx(35.0)
    assert np.isnan(normalized.iloc[2])
    assert np.isnan(normalized.iloc[3])
