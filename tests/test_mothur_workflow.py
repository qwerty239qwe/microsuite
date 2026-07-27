from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.workflows.catalog import WORKFLOWS
from microsuite.workflows.mothur_sop import write_stability_file


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
