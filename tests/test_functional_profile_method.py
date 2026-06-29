from __future__ import annotations

import json
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.functional_profile import functional_profile


def test_tax4fun2_r_script_is_external_asset() -> None:
    packaged_script = files("microsuite.functional.r").joinpath("tax4fun2.R")

    assert packaged_script.is_file()
    text = packaged_script.read_text(encoding="utf-8")
    assert "Tax4Fun2::runRefBlast" in text
    assert "Tax4Fun2::makeFunctionalPrediction" in text


def test_picrust2_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    table = tmp_path / "table.biom"
    rep_seqs = tmp_path / "rep-seqs.fasta"
    table.write_text("placeholder\n", encoding="utf-8")
    rep_seqs.write_text(">s1\nACGT\n", encoding="utf-8")
    calls: list[list[str]] = []

    monkeypatch.setattr(
        shutil,
        "which",
        lambda name: "picrust2_pipeline.py" if name == "picrust2_pipeline.py" else None,
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    functional_profile(
        backend="picrust2",
        table=table,
        rep_seqs=rep_seqs,
        output_dir=tmp_path / "picrust2",
        threads=3,
    )

    assert calls == [
        [
            "picrust2_pipeline.py",
            "-s",
            str(rep_seqs),
            "-i",
            str(table),
            "-o",
            str(tmp_path / "picrust2"),
            "-p",
            "3",
        ]
    ]


def test_tax4fun2_builds_command_and_logs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    table = tmp_path / "otu-table.tsv"
    rep_seqs = tmp_path / "otus.fasta"
    database = tmp_path / "Tax4Fun2_ReferenceData_v2"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    rep_seqs.write_text(">f1\nACGT\n", encoding="utf-8")
    database.mkdir()
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    run_dir = tmp_path / "run"
    functional_profile(
        backend="tax4fun2",
        table=table,
        rep_seqs=rep_seqs,
        database=database,
        output_dir=tmp_path / "tax4fun2",
        threads=2,
        database_mode="Ref100NR",
        min_identity=0.95,
        normalize_pathways=True,
        run_dir=run_dir,
    )

    assert calls
    command = calls[0]
    assert command[0] == "Rscript"
    assert command[1].endswith("microsuite/functional/r/tax4fun2.R") or command[1].endswith(
        "microsuite\\functional\\r\\tax4fun2.R"
    )
    assert command[-4:] == ["2", "Ref100NR", "0.95", "TRUE"]
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "functional_profile"
    assert run["backend"] == "tax4fun2"
    assert run["outputs"] == {"output_dir": str(tmp_path / "tax4fun2")}


def test_humann_builds_command_with_databases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = tmp_path / "reads.fastq.gz"
    nucleotide_db = tmp_path / "chocophlan"
    protein_db = tmp_path / "uniref"
    reads.write_text("placeholder\n", encoding="utf-8")
    nucleotide_db.mkdir()
    protein_db.mkdir()
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "humann" if name == "humann" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    functional_profile(
        backend="humann",
        reads=reads,
        database=nucleotide_db,
        protein_database=protein_db,
        output_dir=tmp_path / "humann",
        threads="4",
    )

    assert calls == [
        [
            "humann",
            "--input",
            str(reads),
            "--output",
            str(tmp_path / "humann"),
            "--threads",
            "4",
            "--nucleotide-database",
            str(nucleotide_db),
            "--protein-database",
            str(protein_db),
        ]
    ]


@pytest.mark.parametrize("backend", ["picrust2", "tax4fun2", "humann"])
def test_functional_profile_reports_missing_external_command(
    backend: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    rep_seqs = tmp_path / "rep-seqs.fasta"
    reads = tmp_path / "reads.fastq"
    database = tmp_path / "db"
    table.write_text("placeholder\n", encoding="utf-8")
    rep_seqs.write_text(">s1\nACGT\n", encoding="utf-8")
    reads.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    database.mkdir()
    monkeypatch.setattr(shutil, "which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError):
        functional_profile(
            backend=backend,
            table=table,
            rep_seqs=rep_seqs,
            reads=reads,
            database=database,
            output_dir=tmp_path / "out",
        )


def test_cli_functional_profile_humann_invokes_method(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "humann" if name == "humann" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    result = CliRunner().invoke(
        app,
        [
            "functional_profile",
            "--backend",
            "humann",
            "--reads",
            str(reads),
            "--output-dir",
            str(tmp_path / "functions"),
            "--threads",
            "2",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert calls[0][:7] == [
        "humann",
        "--input",
        str(reads),
        "--output",
        str(tmp_path / "functions"),
        "--threads",
        "2",
    ]
