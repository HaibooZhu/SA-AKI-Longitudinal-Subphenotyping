from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "09_revision_documents"
    / "audit_reviewer_traceability.py"
)
SPEC = importlib.util.spec_from_file_location("audit_reviewer_traceability", MODULE_PATH)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_response_section_parser_preserves_complete_sections():
    sections = MODULE.parse_sections(
        "# Response\n\n## E.1 First\n\n> comment\n\n**Response:** one\n\n"
        "## E.2 Second\n\n> comment\n\n**Response:** two\n"
    )
    assert set(sections) == {"E.1 First", "E.2 Second"}
    assert "**Response:** one" in sections["E.1 First"]
    assert "**Response:** two" in sections["E.2 Second"]


def test_expected_comment_identifiers_are_unique_and_complete():
    titles = [requirement.title for requirement in MODULE.REQUIREMENTS]
    assert len(titles) == len(set(titles)) == 15
    assert titles[0].startswith("E.1 ")
    assert titles[-1].startswith("R1.m5 ")


def test_token_match_is_case_insensitive_and_fail_closed():
    assert MODULE.contains_all("No External Incremental Value", ("external incremental",))[0]
    passed, missing = MODULE.contains_all("No external value", ("incremental",))
    assert passed is False
    assert missing == ["incremental"]
