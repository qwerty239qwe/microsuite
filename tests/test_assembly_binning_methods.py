from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

import microsuite.api as api
from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.assembly import assemble
from microsuite.methods.binning import bin_contigs


def touch(path: Path, text: str = "placeholder\n") -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_python_sdk_facade_exports_assembly_and_binning() -> None:
    assert api.assemble is assemble
    assert api.bin_contigs is bin_contigs
    assert "assemble" in api.__all__
    assert "bin_contigs" in api.__all__


def test_megahit_builds_paired_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read1 = touch(tmp_path / "R1.fastq.gz")
    read2 = touch(tmp_path / "R2.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "megahit" if name == "megahit" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    assemble(
        backend="megahit",
        read1=read1,
        read2=read2,
        output_dir=tmp_path / "megahit",
        threads=4,
    )

    assert calls == [
        [
            "megahit",
            "-1",
            str(read1),
            "-2",
            str(read2),
            "-o",
            str(tmp_path / "megahit"),
            "-t",
            "4",
        ]
    ]


def test_metaspades_builds_single_command_and_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = touch(tmp_path / "reads.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr(
        shutil, "which", lambda name: "metaspades.py" if name == "metaspades.py" else None
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    run_dir = tmp_path / "run"
    assemble(
        backend="metaspades",
        reads=reads,
        output_dir=tmp_path / "metaspades",
        threads=2,
        run_dir=run_dir,
    )

    assert calls == [
        [
            "metaspades.py",
            "-s",
            str(reads),
            "-o",
            str(tmp_path / "metaspades"),
            "-t",
            "2",
        ]
    ]
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "assemble"
    assert run["backend"] == "metaspades"
    assert run["outputs"] == {"output_dir": str(tmp_path / "metaspades")}


def test_idba_ud_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reads = touch(tmp_path / "reads.fa", ">r1\nACGT\n")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "idba_ud" if name == "idba_ud" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    assemble(
        backend="idba-ud",
        reads=reads,
        output_dir=tmp_path / "idba",
        threads=3,
    )

    assert calls == [
        [
            "idba_ud",
            "-r",
            str(reads),
            "-o",
            str(tmp_path / "idba"),
            "--num_threads",
            "3",
        ]
    ]


def test_assembly_validates_inputs_and_missing_tool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read1 = touch(tmp_path / "R1.fastq.gz")
    reads = touch(tmp_path / "reads.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="either --reads"):
        assemble(
            backend="megahit",
            read1=read1,
            reads=reads,
            output_dir=tmp_path / "out",
        )

    monkeypatch.setattr(shutil, "which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="MEGAHIT"):
        assemble(backend="megahit", read1=read1, output_dir=tmp_path / "out")


def test_metabat2_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contigs = touch(tmp_path / "contigs.fa", ">c1\nACGT\n")
    depth = touch(tmp_path / "depth.tsv")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "metabat2" if name == "metabat2" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    bin_contigs(
        backend="metabat2",
        contigs=contigs,
        depth=depth,
        output_dir=tmp_path / "metabat",
        prefix="sample-bin",
        threads=5,
    )

    assert calls == [
        [
            "metabat2",
            "-i",
            str(contigs),
            "-a",
            str(depth),
            "-o",
            str(tmp_path / "metabat" / "sample-bin"),
            "-t",
            "5",
        ]
    ]


def test_maxbin2_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contigs = touch(tmp_path / "contigs.fa", ">c1\nACGT\n")
    abundance = touch(tmp_path / "abundance.tsv")
    calls: list[list[str]] = []

    monkeypatch.setattr(
        shutil, "which", lambda name: "run_MaxBin.pl" if name == "run_MaxBin.pl" else None
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    bin_contigs(
        backend="maxbin2",
        contigs=contigs,
        abundance=abundance,
        output_dir=tmp_path / "maxbin",
        threads=6,
    )

    assert calls == [
        [
            "run_MaxBin.pl",
            "-contig",
            str(contigs),
            "-abund",
            str(abundance),
            "-out",
            str(tmp_path / "maxbin" / "bin"),
            "-thread",
            "6",
        ]
    ]


def test_concoct_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    contigs = touch(tmp_path / "contigs.fa", ">c1\nACGT\n")
    coverage = touch(tmp_path / "coverage.tsv")
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "concoct" if name == "concoct" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    bin_contigs(
        backend="concoct",
        contigs=contigs,
        coverage=coverage,
        output_dir=tmp_path / "concoct",
        threads=7,
    )

    assert calls == [
        [
            "concoct",
            "--composition_file",
            str(contigs),
            "--coverage_file",
            str(coverage),
            "-b",
            str(tmp_path / "concoct"),
            "-t",
            "7",
        ]
    ]


@pytest.mark.parametrize(
    ("backend", "missing"),
    [("metabat2", "--depth"), ("maxbin2", "--abundance"), ("concoct", "--coverage")],
)
def test_binning_requires_backend_specific_abundance_inputs(
    backend: str, missing: str, tmp_path: Path
) -> None:
    contigs = touch(tmp_path / "contigs.fa", ">c1\nACGT\n")

    with pytest.raises(MicrobiomeSuiteError, match=missing):
        bin_contigs(backend=backend, contigs=contigs, output_dir=tmp_path / "bins")


def test_cli_assemble_and_bin_help_and_invocation(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = touch(tmp_path / "reads.fastq.gz")
    contigs = touch(tmp_path / "contigs.fa", ">c1\nACGT\n")
    depth = touch(tmp_path / "depth.tsv")
    calls: list[list[str]] = []

    def fake_which(name: str) -> str | None:
        return {"megahit": "megahit", "metabat2": "metabat2"}.get(name)

    monkeypatch.setattr(shutil, "which", fake_which)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    runner = CliRunner()
    methods = runner.invoke(app, ["methods"])
    assert methods.exit_code == 0
    assert "assemble" in methods.stdout
    assert "megahit" in methods.stdout
    assert "bin" in methods.stdout
    assert "metabat2" in methods.stdout

    assemble_help = runner.invoke(app, ["assemble", "--help"])
    assert assemble_help.exit_code == 0
    assert "--read1" in assemble_help.stdout
    assert "--reads" in assemble_help.stdout

    bin_help = runner.invoke(app, ["bin", "--help"])
    assert bin_help.exit_code == 0
    assert "--contigs" in bin_help.stdout
    assert "--depth" in bin_help.stdout

    result = runner.invoke(
        app,
        [
            "assemble",
            "--backend",
            "megahit",
            "--reads",
            str(reads),
            "--output-dir",
            str(tmp_path / "assembly"),
        ],
    )
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(
        app,
        [
            "bin",
            "--backend",
            "metabat2",
            "--contigs",
            str(contigs),
            "--depth",
            str(depth),
            "--output-dir",
            str(tmp_path / "bins"),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert calls[0][:3] == ["megahit", "-r", str(reads)]
    assert calls[1][:5] == ["metabat2", "-i", str(contigs), "-a", str(depth)]
