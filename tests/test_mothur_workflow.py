from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.workflows.catalog import WORKFLOWS
from microsuite.workflows.mothur_sop import run_mothur_sop, write_stability_file


def test_mothur_workflow_is_in_the_catalog() -> None:
    assert any(workflow.name == "mothur" for workflow in WORKFLOWS)


def test_write_stability_file_pairs_r1_and_r2(tmp_path: Path) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    for name in (
        "sampleA_R1.fastq.gz",
        "sampleA_R2.fastq.gz",
        "sampleB_R1.fastq.gz",
        "sampleB_R2.fastq.gz",
    ):
        (reads / name).write_text("", encoding="utf-8")

    stability = write_stability_file(reads, tmp_path / "stability.files")

    lines = [line.split("\t") for line in stability.read_text(encoding="utf-8").splitlines()]
    assert [line[0] for line in lines] == ["sampleA", "sampleB"]
    assert lines[0][1].endswith("sampleA_R1.fastq.gz")
    assert lines[0][2].endswith("sampleA_R2.fastq.gz")


def test_write_stability_file_rejects_unpaired_reads(tmp_path: Path) -> None:
    # A dropped mate silently halves a sample's depth if it is not caught here.
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sampleA_R1.fastq.gz").write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="sampleA"):
        write_stability_file(reads, tmp_path / "stability.files")


def test_write_stability_file_rejects_empty_directory(tmp_path: Path) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()

    with pytest.raises(MicrobiomeSuiteError, match="No paired FASTQ"):
        write_stability_file(reads, tmp_path / "stability.files")


@pytest.mark.parametrize(
    ("r1_name", "r2_name", "expected_sample"),
    [
        ("sampleA_R1.fastq.gz", "sampleA_R2.fastq.gz", "sampleA"),
        (
            "sampleA_S1_L001_R1_001.fastq.gz",
            "sampleA_S1_L001_R2_001.fastq.gz",
            "sampleA",
        ),
        ("sampleA_r1.fastq.gz", "sampleA_r2.fastq.gz", "sampleA"),
        ("sampleA.R1.fastq.gz", "sampleA.R2.fastq.gz", "sampleA"),
        ("sample_A_R1_001.fastq.gz", "sample_A_R2_001.fastq.gz", "sample_A"),
    ],
    ids=[
        "plain_r1_r2",
        "illumina_bcl2fastq_lane_naming",
        "lowercase_mate_marker",
        "dot_separated_mate",
        "underscore_in_sample_name",
    ],
)
def test_write_stability_file_parses_known_filename_formats(
    tmp_path: Path, r1_name: str, r2_name: str, expected_sample: str
) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / r1_name).write_text("", encoding="utf-8")
    (reads / r2_name).write_text("", encoding="utf-8")

    stability = write_stability_file(reads, tmp_path / "stability.files")
    lines = [line.split("\t") for line in stability.read_text(encoding="utf-8").splitlines()]

    assert len(lines) == 1, f"{r1_name}/{r2_name}: expected exactly one sample"
    assert lines[0][0] == expected_sample, f"{r1_name}/{r2_name}: wrong sample name"
    assert lines[0][1].endswith(r1_name), f"{r1_name}/{r2_name}: wrong R1 path"
    assert lines[0][2].endswith(r2_name), f"{r1_name}/{r2_name}: wrong R2 path"


def test_write_stability_file_rejects_multiple_lanes_for_one_sample(tmp_path: Path) -> None:
    # bcl2fastq splits one sample across lanes; mothur's stability file can only
    # hold one R1/R2 pair per sample, so this must raise instead of silently
    # dropping a lane's reads.
    reads = tmp_path / "reads"
    reads.mkdir()
    for name in (
        "sampleA_S1_L001_R1_001.fastq.gz",
        "sampleA_S1_L001_R2_001.fastq.gz",
        "sampleA_S1_L002_R1_001.fastq.gz",
        "sampleA_S1_L002_R2_001.fastq.gz",
    ):
        (reads / name).write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="sampleA"):
        write_stability_file(reads, tmp_path / "stability.files")


def test_write_stability_file_ignores_non_fastq_files(tmp_path: Path) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    for name in (
        "sampleA_R1.fastq.gz",
        "sampleA_R2.fastq.gz",
        "notes.txt",
        "README.md",
    ):
        (reads / name).write_text("", encoding="utf-8")

    stability = write_stability_file(reads, tmp_path / "stability.files")
    lines = stability.read_text(encoding="utf-8").splitlines()

    assert len(lines) == 1
    assert lines[0].startswith("sampleA\t")


def test_write_stability_file_rejects_unmatched_fastq_file(tmp_path: Path) -> None:
    # A FASTQ file the mate regex cannot place is a sample about to go missing;
    # it must raise, not be skipped like a genuinely non-FASTQ file.
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sampleA_R1.fastq.gz").write_text("", encoding="utf-8")
    (reads / "sampleA_R2.fastq.gz").write_text("", encoding="utf-8")
    (reads / "weirdfile.fastq.gz").write_text("", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="weirdfile.fastq.gz"):
        write_stability_file(reads, tmp_path / "stability.files")


def test_run_mothur_sop_threads_files_correctly_between_steps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Every existing test targets write_stability_file or catalog membership;
    # nothing pinned the hand-offs between make.contigs -> cluster ->
    # tax_classify. An earlier task in this feature shipped a mis-threaded
    # file with a fully green suite. This monkeypatches run_mothur, cluster,
    # and tax_classify at the module boundary and asserts the exact paths
    # crossing each hand-off.
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sampleA_R1.fastq.gz").write_text("", encoding="utf-8")
    (reads / "sampleA_R2.fastq.gz").write_text("", encoding="utf-8")

    output_dir = tmp_path / "out"
    reference_alignment = tmp_path / "silva.align"
    reference_alignment.write_text(">ref\nAC-GT\n", encoding="utf-8")
    taxonomy_reference = tmp_path / "trainset.fasta"
    taxonomy_reference.write_text(">r\nACGT\n", encoding="utf-8")
    taxonomy_map = tmp_path / "trainset.tax"
    taxonomy_map.write_text("r\tBacteria;\n", encoding="utf-8")

    contigs_fasta = output_dir / "contigs" / "stability.trim.contigs.fasta"
    scrap_fasta = output_dir / "contigs" / "stability.scrap.contigs.fasta"
    contigs_count_table = output_dir / "contigs" / "stability.contigs.count_table"

    run_mothur_calls: list[dict[str, object]] = []
    cluster_calls: list[dict[str, object]] = []
    tax_classify_calls: list[dict[str, object]] = []

    def fake_run_mothur(command: str, params: dict[str, str], **kwargs: object) -> list[Path]:
        run_mothur_calls.append({"command": command, "params": params, **kwargs})
        # make.contigs emits the assembled FASTA, the scrap FASTA, and the
        # count table (the ONLY carrier of read->sample identity) in one block.
        return [contigs_fasta, scrap_fasta, contigs_count_table]

    def fake_cluster(**kwargs: object) -> None:
        cluster_calls.append(kwargs)

    def fake_tax_classify(**kwargs: object) -> None:
        tax_classify_calls.append(kwargs)

    monkeypatch.setattr("microsuite.workflows.mothur_sop.run_mothur", fake_run_mothur)
    monkeypatch.setattr("microsuite.workflows.mothur_sop.cluster", fake_cluster)
    monkeypatch.setattr("microsuite.workflows.mothur_sop.tax_classify", fake_tax_classify)

    run_mothur_sop(
        reads_dir=reads,
        output_dir=output_dir,
        reference_alignment=reference_alignment,
        taxonomy_reference=taxonomy_reference,
        taxonomy_map=taxonomy_map,
    )

    assert len(cluster_calls) == 1
    assert len(tax_classify_calls) == 1

    # 1. The contigs FASTA from make.contigs (NOT the scrap FASTA) reaches
    # cluster(rep_seqs=...).
    assert cluster_calls[0]["rep_seqs"] == contigs_fasta

    # 2. cluster's output_rep_seqs becomes tax_classify's rep_seqs.
    assert tax_classify_calls[0]["rep_seqs"] == cluster_calls[0]["output_rep_seqs"]

    # 3. cluster's output_otu_list/output_count_table become tax_classify's
    # otu_list/count_table -- the Finding 1 per-OTU consensus wiring.
    assert tax_classify_calls[0]["otu_list"] == cluster_calls[0]["output_otu_list"]
    assert tax_classify_calls[0]["count_table"] == cluster_calls[0]["output_count_table"]
    assert cluster_calls[0]["output_otu_list"] is not None
    assert cluster_calls[0]["output_count_table"] is not None

    # 4. make.contigs's .count_table (NOT the trim or scrap FASTA) reaches
    # cluster(count_table=...). It is the ONLY carrier of read->sample
    # identity in the whole pipeline -- mothur never renames reads to encode
    # a sample, so dropping this makes every downstream table single-sample
    # regardless of how many samples were fed in.
    assert cluster_calls[0]["count_table"] == contigs_count_table
