from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "03_cluster_robustness"
    / "audit_mixak_provenance.py"
)
SPEC = importlib.util.spec_from_file_location("audit_mixak_provenance", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_notebook_scanner_reads_code_cells_only(tmp_path):
    notebook = {
        "cells": [
            {"cell_type": "code", "source": ["prior.b=list(Kmax=3)\n"]},
            {"cell_type": "markdown", "source": ["Kmax=8\n"]},
            {"cell_type": "code", "source": ["x <- quant.comp.prob[[\"50%\"]] / 2\n"]},
        ]
    }
    path = tmp_path / "model.ipynb"
    path.write_text(json.dumps(notebook), encoding="utf-8")
    text = MODULE.source_text(path)
    assert MODULE.K_PATTERN.findall(text) == ["3"]
    assert len(MODULE.PROBABILITY_DIVISION_PATTERN.findall(text)) == 1


def test_python_scanner_finds_candidate_k_and_probability_scaling(tmp_path):
    path = tmp_path / "model.py"
    path.write_text(
        'code = "prior.b = list(Kmax = 5)"\n'
        'other = "mod$quant.comp.prob[[\\\"50%\\\"]] / 2"\n',
        encoding="utf-8",
    )
    inventory = MODULE.audit_code([tmp_path])
    assert json.loads(inventory.loc[0, "candidate_k_literals"]) == [5]
    assert inventory.loc[0, "posterior_probability_division_by_2_hits"] == 1


def test_python_source_does_not_claim_opaque_rdata_objects_are_verified():
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert '"verified_model_objects"' not in source
    assert '"expected_model_objects"' in source


def test_provenance_status_is_not_left_as_an_unfinished_table_action():
    source = MODULE_PATH.read_text(encoding="utf-8")
    assert "PASS_WITH_REQUIRED_TABLE_CORRECTION" not in source
    assert "PASS_TRACEABLE_K2_K5_ONLY" in source
