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
    assert "Tax4Fun2::makeFunctionalPrediction" in text
    assert 'Sys.which("makeblastdb")' in text
    assert "tax4fun2_manifest.json" in text


def _write_tax4fun2_database(path: Path, mode: str = "Ref99NR") -> None:
    profiles = path / mode
    kegg = path / "KEGG"
    profiles.mkdir(parents=True)
    kegg.mkdir(parents=True)
    (profiles / f"{mode}.fasta").write_text(">REF\nACGT\n", encoding="utf-8")
    (profiles / "REF.tbl.gz").write_bytes(b"placeholder")
    (kegg / "ko.txt").write_text("ko\tdescription\tptw_count\n", encoding="utf-8")
    (kegg / "ko2ptw.txt").write_text("nrow\tptw\n", encoding="utf-8")
    (kegg / "ptw.txt").write_text("ptw\tdescription\n", encoding="utf-8")


def _write_tax4fun2_outputs(path: Path) -> None:
    path.mkdir(parents=True)
    (path / "functional_prediction.tsv").write_text(
        "KO\ts1\tdescription\nK00001\t1\ttest\n", encoding="utf-8"
    )
    (path / "pathway_prediction.tsv").write_text(
        "pathway\ts1\tdescription\nmap00010\t1\ttest\n", encoding="utf-8"
    )
    (path / "coverage.tsv").write_text(
        "sample\tfeature_fraction_used\tsequence_fraction_used\ns1\t1\t1\n",
        encoding="utf-8",
    )
    (path / "tax4fun2_manifest.json").write_text(
        '{"schema_version":"microsuite-tax4fun2.v1"}\n', encoding="utf-8"
    )


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
    _write_tax4fun2_database(database, "Ref100NR")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        _write_tax4fun2_outputs(Path(command[-1]))
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)

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
    assert command[-5:-1] == ["2", "Ref100NR", "0.95", "TRUE"]
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "functional_profile"
    assert run["backend"] == "tax4fun2"
    assert run["outputs"] == {
        "coverage": str(tmp_path / "tax4fun2" / "coverage.tsv"),
        "functions": str(tmp_path / "tax4fun2" / "functional_prediction.tsv"),
        "manifest": str(tmp_path / "tax4fun2" / "tax4fun2_manifest.json"),
        "pathways": str(tmp_path / "tax4fun2" / "pathway_prediction.tsv"),
    }
    assert run["params"]["tax4fun2_version"] == "1.1.5"
    assert (tmp_path / "tax4fun2" / "coverage.tsv").is_file()


@pytest.mark.parametrize(
    ("table_text", "fasta_text", "message"),
    [
        ("feature-id\ts1\nf1\t-1\n", ">f1\nACGT\n", "finite and non-negative"),
        ("feature-id\ts1\nf1\tbad\n", ">f1\nACGT\n", "not numeric"),
        ("feature-id\ts1\nf1\t0\n", ">f1\nACGT\n", "positive total abundance"),
        ("feature-id\ts1\nf1\t1\n", ">other\nACGT\n", "must match exactly"),
        ("feature-id\ts1\nf1\t1\nf1\t2\n", ">f1\nACGT\n", "duplicated"),
    ],
)
def test_tax4fun2_rejects_invalid_inputs(
    tmp_path: Path, table_text: str, fasta_text: str, message: str
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text(table_text, encoding="utf-8")
    fasta.write_text(fasta_text, encoding="utf-8")
    _write_tax4fun2_database(database)

    with pytest.raises(MicrobiomeSuiteError, match=message):
        functional_profile(
            backend="tax4fun2",
            table=table,
            rep_seqs=fasta,
            database=database,
            output_dir=tmp_path / "out",
        )


def test_tax4fun2_rejects_incomplete_reference(tmp_path: Path) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    database.mkdir()

    with pytest.raises(MicrobiomeSuiteError, match="reference data is incomplete"):
        functional_profile(
            backend="tax4fun2",
            table=table,
            rep_seqs=fasta,
            database=database,
            output_dir=tmp_path / "out",
        )


def test_tax4fun2_builds_docker_command_and_keeps_stable_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database)
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: f"/usr/bin/{name}")

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        if "run" in command:
            output_arg = command[-1]
            staged = next(tmp_path.parent.glob(".microsuite-tax4fun2-*/result"), None)
            if staged is None:
                stage_root = next(tmp_path.glob(".microsuite-tax4fun2-*"))
                staged = stage_root / Path(output_arg).name
            _write_tax4fun2_outputs(staged)
            return subprocess.CompletedProcess(command, 0, "ok\n", "")
        return subprocess.CompletedProcess(command, 0, "sha256:test\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    output = tmp_path / "out"
    functional_profile(
        backend="tax4fun2",
        table=table,
        rep_seqs=fasta,
        database=database,
        output_dir=output,
        runtime="docker",
        image="example/tax4fun2:1.1.5",
    )

    docker = calls[0]
    assert docker[:3] == ["docker", "run", "--rm"]
    assert "example/tax4fun2:1.1.5" in docker
    assert "/opt/microsuite/tax4fun2.R" in docker
    assert (output / "functional_prediction.tsv").is_file()
    assert (output / "tax4fun2_container.json").is_file()


def test_tax4fun2_failed_run_does_not_replace_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    output = tmp_path / "out"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database)
    output.mkdir()
    marker = output / "keep.txt"
    marker.write_text("original\n", encoding="utf-8")
    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 1, "", "failed"),
    )

    with pytest.raises(MicrobiomeSuiteError, match="failed"):
        functional_profile(
            backend="tax4fun2",
            table=table,
            rep_seqs=fasta,
            database=database,
            output_dir=output,
            force=True,
        )

    assert marker.read_text(encoding="utf-8") == "original\n"


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
    table.write_text("feature-id\ts1\ns1\t1\n", encoding="utf-8")
    rep_seqs.write_text(">s1\nACGT\n", encoding="utf-8")
    reads.write_text("@r1\nACGT\n+\n!!!!\n", encoding="utf-8")
    if backend == "tax4fun2":
        _write_tax4fun2_database(database)
    else:
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


def test_cli_functional_profile_tax4fun2_invokes_hardened_backend(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.tsv"
    fasta = tmp_path / "rep-seqs.fasta"
    database = tmp_path / "reference"
    table.write_text("feature-id\ts1\nf1\t1\n", encoding="utf-8")
    fasta.write_text(">f1\nACGT\n", encoding="utf-8")
    _write_tax4fun2_database(database)
    calls: list[list[str]] = []
    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        _write_tax4fun2_outputs(Path(command[-1]))
        return subprocess.CompletedProcess(command, 0, "ok\n", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    result = CliRunner().invoke(
        app,
        [
            "functional_profile",
            "--backend",
            "tax4fun2",
            "--table",
            str(table),
            "--rep-seqs",
            str(fasta),
            "--database",
            str(database),
            "--output-dir",
            str(tmp_path / "functions"),
            "--min-identity",
            "0.95",
            "--normalize-pathways",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert calls[0][-5:-1] == ["1", "Ref99NR", "0.95", "TRUE"]
    assert (tmp_path / "functions" / "functional_prediction.tsv").is_file()
