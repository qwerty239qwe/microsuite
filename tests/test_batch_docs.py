from __future__ import annotations

from pathlib import Path

from microsuite.batch.backends import BATCH_BACKENDS

ROOT = Path(__file__).resolve().parents[1]


def test_every_backend_is_documented() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    for name in BATCH_BACKENDS:
        assert name in text, f"{name} is not documented"


def test_the_leakage_hazard_is_stated() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    assert "--target-col" in text
    assert "plsda-batch" in text
    # The document must say what goes wrong, not merely that an option exists.
    assert "inflat" in text.lower()


def test_the_scale_contract_is_documented() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    for value_type in ("counts", "relative", "clr"):
        assert value_type in text


def test_correction_is_not_presented_as_a_substitute_for_modelling() -> None:
    text = (ROOT / "docs" / "batch_correction.md").read_text(encoding="utf-8")
    assert "covariate" in text.lower()
