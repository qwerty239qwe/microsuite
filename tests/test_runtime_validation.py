from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.validation import validate_output_file, validate_outputs


def test_missing_output_raises(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="not created"):
        validate_output_file(tmp_path / "nope.tsv")


def test_empty_output_raises_with_placeholder_hint(tmp_path: Path) -> None:
    p = tmp_path / "empty.tsv"
    p.write_bytes(b"")
    with pytest.raises(MicrobiomeSuiteError, match="empty"):
        validate_output_file(p)
    # message mentions the cloud-placeholder possibility
    try:
        validate_output_file(p)
    except MicrobiomeSuiteError as exc:
        assert "placeholder" in str(exc).lower()


def test_nonempty_tsv_passes(tmp_path: Path) -> None:
    p = tmp_path / "t.tsv"
    p.write_text("a\tb\n1\t2\n", encoding="utf-8")
    validate_output_file(p)  # no raise


def test_valid_gzip_passes(tmp_path: Path) -> None:
    p = tmp_path / "reads.fastq.gz"
    p.write_bytes(gzip.compress(b"@r\nACGT\n+\nIIII\n"))
    validate_output_file(p)  # no raise


def test_invalid_gzip_raises(tmp_path: Path) -> None:
    p = tmp_path / "bad.fastq.gz"
    p.write_bytes(b"this is not gzip but has content")
    with pytest.raises(MicrobiomeSuiteError, match="gzip"):
        validate_output_file(p)


def test_allow_empty(tmp_path: Path) -> None:
    p = tmp_path / "empty.txt"
    p.write_bytes(b"")
    validate_output_file(p, allow_empty=True)  # no raise


def test_validate_outputs_raises_on_first_bad(tmp_path: Path) -> None:
    good = tmp_path / "g.tsv"
    good.write_text("x\n", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        validate_outputs({"good": str(good), "missing": str(tmp_path / "no.tsv")})
