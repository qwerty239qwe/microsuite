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


# Pipeline-configuration fixtures (captured 2026-07-28): a real 2-sample
# paired-end run through mothur 1.48.5, aligned and with a count table
# carrying sample groups. See tests/fixtures/mothur/README.md -- four mothur
# commands behave differently here than on the unaligned, group-less toy
# input the fixtures above were captured on.


def test_select_output_picks_make_contigs_count_table_alongside_scrap_fasta() -> None:
    # make.contigs's .count_table carries group columns and is the ONLY
    # carrier of read->sample identity; make.contigs does not rename reads.
    outputs = parse_mothur_outputs(_fixture("make_contigs_paired.txt"))

    fasta = select_output(outputs, ".fasta", step="make.contigs", exclude=("scrap",))
    count_table = select_output(outputs, ".count_table", step="make.contigs")

    assert fasta.name == "stability.trim.contigs.fasta"
    assert count_table.name == "stability.contigs.count_table"


def test_select_output_picks_screen_seqs_align_when_nothing_removed() -> None:
    # Input here is align.seqs's .align, so screen.seqs's output is
    # .good.align, never .good.fasta. Nothing was screened out, so no
    # .count_table is emitted at all.
    outputs = parse_mothur_outputs(_fixture("screen_seqs_aligned.txt"))

    chosen = select_output(outputs, ".align", step="screen.seqs")

    assert chosen.name == "stability.trim.contigs.unique.good.align"
    assert not any(path.name.endswith(".count_table") for path in outputs)


def test_select_output_picks_screen_seqs_count_table_when_removed() -> None:
    # When sequences ARE removed, screen.seqs emits .good.align, .bad.accnos,
    # AND .good.count_table. Two "Output File Names:" blocks appear (an
    # internal remove.seqs runs first) -- last-block-wins keeps this
    # unambiguous.
    outputs = parse_mothur_outputs(_fixture("screen_seqs_removed.txt"))

    aligned = select_output(outputs, ".align", step="screen.seqs")
    count_table = select_output(outputs, ".count_table", step="screen.seqs")
    accnos = select_output(outputs, ".accnos", step="screen.seqs")

    assert aligned.name == "q.unique.good.align"
    assert count_table.name == "q.good.count_table"
    assert accnos.name == "q.unique.bad.accnos"
    assert "q.pick.count_table" not in [p.name for p in outputs]


def test_select_output_chimera_vsearch_grouped_has_no_fasta() -> None:
    # With a GROUPED count table, chimera.vsearch emits .count_table +
    # .accnos + .chimeras and NO .fasta -- the inverse of chimera_vsearch.txt
    # (group-less input), which emits .fasta and no .count_table. Feeding a
    # caller the wrong one of these is exactly the D3 defect.
    outputs = parse_mothur_outputs(_fixture("chimera_vsearch_grouped.txt"))

    count_table = select_output(outputs, ".count_table", step="chimera.vsearch")
    accnos = select_output(outputs, ".accnos", step="chimera.vsearch")

    assert count_table.name == (
        "stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.count_table"
    )
    assert accnos.name == (
        "stability.trim.contigs.unique.good.filter.unique.precluster.denovo.vsearch.accnos"
    )
    with pytest.raises(MicrobiomeSuiteError, match="no '.fasta' output"):
        select_output(outputs, ".fasta", step="chimera.vsearch")


def test_select_output_cluster_picks_opti_mcc_list() -> None:
    # Real cluster output is named opti_mcc, not opti_tptn -- opti_tptn only
    # appears when the distance file is blank and cluster aborts.
    outputs = parse_mothur_outputs(_fixture("cluster_opti.txt"))

    chosen = select_output(outputs, ".list", step="cluster")

    assert "opti_mcc" in chosen.name
    assert "opti_tptn" not in chosen.name


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
