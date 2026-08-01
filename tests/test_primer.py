from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from typer.testing import CliRunner

from microsuite.cli.app import app
from microsuite.primer import check_fastq_primers, primer_check_fails


def write_fastq(path: Path, sequences: list[str]) -> None:
    text = "".join(
        f"@read-{i}\n{sequence}\n+\n{'I' * len(sequence)}\n" for i, sequence in enumerate(sequences)
    )
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        path.write_text(text, encoding="utf-8")


def test_primer_checker_handles_iupac_and_key_barcode(tmp_path: Path) -> None:
    read = tmp_path / "reads.fastq.gz"
    write_fastq(
        read,
        [
            "TCAGACGTACGTACGGATTAGATACCCTGGTAGTACGT",
            "TCAGTTTTCCCCAAGGATTAGATACCCTGGTAGTACGT",
        ],
    )

    report = check_fastq_primers(
        [("R1", read)],
        cutadapt={"front": "^TCAGNNNNNNNNNNGGATTAGATACCCTGGTAGT"},
        primer_check={"max_mismatches": 0, "min_match_rate": 1.0},
    )

    assert report["status"] == "passed"
    assert report["pattern_results"]["front"]["match_rate"] == 1.0
    assert report["records_examined"] == 2


def test_primer_checker_reports_low_match_rate_without_failing_warn_mode(tmp_path: Path) -> None:
    read = tmp_path / "reads.fastq"
    write_fastq(read, ["ACGTACGTACGT"])

    report = check_fastq_primers(
        [("single", read)],
        cutadapt={"adapter": "CTGAGCCAGGATCAAACTCT"},
        primer_check={"mode": "warn", "min_match_rate": 0.8},
    )

    assert report["status"] == "warning"
    assert primer_check_fails(report, "warn") is False
    assert primer_check_fails(report, "error") is True


def test_primer_checker_writes_json_from_cli(tmp_path: Path) -> None:
    read = tmp_path / "reads.fastq"
    output = tmp_path / "primer-check.json"
    write_fastq(read, ["ACGTCTGAGCCAGGATCAAACTCT"])

    result = CliRunner().invoke(
        app,
        [
            "primer-check",
            "--input",
            str(read),
            "--adapter",
            "CTGAGCCAGGATCAAACTCT",
            "--min-match-rate",
            "1",
            "--output",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert json.loads(output.read_text(encoding="utf-8"))["status"] == "passed"


def test_primer_checker_rejects_invalid_iupac() -> None:
    with pytest.raises(ValueError, match="unsupported bases"):
        check_fastq_primers(
            [],
            cutadapt={"front": "^ACGTZ"},
            primer_check={},
        )
