from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.mothur import (
    check_mothur_errors,
    parse_mothur_outputs,
    select_output,
)

FIXTURES = Path(__file__).resolve().parent / "fixtures" / "mothur"


def _fixture(name: str) -> str:
    return (FIXTURES / name).read_text(encoding="utf-8")


def test_parse_returns_output_paths_in_order() -> None:
    outputs = parse_mothur_outputs(_fixture("unique_seqs.txt"))

    assert [p.name for p in outputs] == [
        "test.unique.fasta",
        "test.count_table",
    ]


def test_parse_returns_empty_when_block_absent() -> None:
    assert parse_mothur_outputs("mothur > quit()\n") == []


def test_select_output_matches_by_suffix() -> None:
    outputs = parse_mothur_outputs(_fixture("unique_seqs.txt"))

    assert select_output(outputs, ".fasta", step="unique.seqs").name.endswith(".unique.fasta")
    assert select_output(outputs, ".count_table", step="unique.seqs").name.endswith(".count_table")


def test_select_output_raises_when_missing_and_lists_what_was_produced() -> None:
    outputs = parse_mothur_outputs(_fixture("unique_seqs.txt"))

    with pytest.raises(MicrobiomeSuiteError) as excinfo:
        select_output(outputs, ".shared", step="unique.seqs")

    message = str(excinfo.value)
    assert "unique.seqs" in message
    assert ".shared" in message
    # The message must name what mothur actually produced, or the user is blind.
    assert "test.unique.fasta" in message


def test_select_output_raises_on_ambiguous_match() -> None:
    # make.contigs emits both .trim.contigs.fasta and .scrap.contigs.fasta.
    # Silently taking the first would hand the scrap reads downstream.
    outputs = parse_mothur_outputs(_fixture("make_contigs.txt"))

    with pytest.raises(MicrobiomeSuiteError, match="ambiguous"):
        select_output(outputs, ".fasta", step="make.contigs")


def test_select_output_exclude_resolves_ambiguity() -> None:
    outputs = parse_mothur_outputs(_fixture("make_contigs.txt"))

    chosen = select_output(outputs, ".fasta", step="make.contigs", exclude=("scrap",))

    assert chosen.name == "stability.trim.contigs.fasta"


def test_select_output_raises_when_no_outputs_at_all() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="produced no output files"):
        select_output([], ".fasta", step="screen.seqs")


def test_check_errors_raises_on_anchored_error_line() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="did not complete align.seqs"):
        check_mothur_errors(_fixture("error_on_failure.txt"), step="align.seqs")


def test_check_errors_ignores_the_summary_banner() -> None:
    # mothur closes a failed run with "Detected 1 [ERROR] messages, please review."
    # A bare '[ERROR]' substring scan matches that banner too and reports one
    # failure twice. The anchor is "[ERROR]: " with the colon and space.
    message = str(
        pytest.raises(
            MicrobiomeSuiteError,
            check_mothur_errors,
            _fixture("error_on_failure.txt"),
            step="align.seqs",
        ).value
    )

    assert "Detected" not in message
    assert message.count("[ERROR]") == 1


def test_check_errors_passes_clean_output() -> None:
    check_mothur_errors(_fixture("unique_seqs.txt"), step="unique.seqs")


def test_check_errors_ignores_warnings() -> None:
    check_mothur_errors("[WARNING]: blank sequence removed\n", step="screen.seqs")
