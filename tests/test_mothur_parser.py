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


def test_parse_takes_last_block_when_script_runs_multiple_commands() -> None:
    # multi_block.txt captures `unique.seqs` followed by `summary.seqs` in one
    # script. unique.seqs emits its own "Output File Names:" block (test.unique.fasta,
    # test.count_table) before summary.seqs runs and emits a second one
    # (test.unique.summary). If parse_mothur_outputs stopped at the first block, the
    # caller asking for summary.seqs's output would silently get unique.seqs's files
    # instead — the preamble command's outputs must not shadow the command the
    # caller actually asked for. Last-block-wins is what makes that safe.
    outputs = parse_mothur_outputs(_fixture("multi_block.txt"))

    assert [p.name for p in outputs] == ["test.unique.summary"]


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


def test_select_output_picks_taxonomy_not_tax_summary() -> None:
    # classify.seqs emits both a .wang.taxonomy and a .wang.tax.summary file.
    # The near-miss is the entire point of pinning this fixture: a caller
    # asking for ".taxonomy" must not be handed the summary file instead.
    outputs = parse_mothur_outputs(_fixture("classify_seqs.txt"))

    chosen = select_output(outputs, ".taxonomy", step="classify.seqs")

    assert chosen.name == "q.unique.ref.wang.taxonomy"
    assert "tax.summary" not in chosen.name


def test_select_output_picks_oturep_fasta_alongside_count_table() -> None:
    # get.oturep (run without the invalid `label` param) emits a
    # .rep.count_table and a .rep.fasta in the same block; select_output must
    # resolve ".fasta" to the fasta file without the count table interfering.
    outputs = parse_mothur_outputs(_fixture("get_oturep.txt"))

    chosen = select_output(outputs, ".fasta", step="get.oturep")

    assert chosen.name == "q.unique.opti_tptn.0.03.0.03.rep.fasta"
    assert "q.unique.opti_tptn.0.03.0.03.rep.count_table" in [p.name for p in outputs]


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
