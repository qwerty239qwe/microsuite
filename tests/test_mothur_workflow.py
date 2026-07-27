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
