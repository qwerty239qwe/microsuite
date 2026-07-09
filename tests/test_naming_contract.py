from __future__ import annotations

from pathlib import Path

import pytest
from naming_contract_cases import CASES, NamingCase

from microsuite.methods.denoise import _expected_sample_ids


@pytest.mark.parametrize("case", CASES, ids=lambda c: c.label)
def test_expected_sample_ids_matches_contract(case: NamingCase, tmp_path: Path) -> None:
    for name in case.filenames:
        (tmp_path / name).write_text("x", encoding="utf-8")
    result = _expected_sample_ids(tmp_path, paired=case.paired)
    assert result == set(case.expected)
