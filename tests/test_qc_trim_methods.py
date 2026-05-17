from __future__ import annotations

import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.qc import qc
from microsuite.methods.trim import trim


def touch(path: Path) -> Path:
    path.write_text("placeholder", encoding="utf-8")
    return path


def test_qc_fastqc_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "fastqc" if name == "fastqc" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc(backend="fastqc", inputs=[read], output_dir=tmp_path / "qc", threads=2)

    assert calls == [["fastqc", "--outdir", str(tmp_path / "qc"), "--threads", "2", str(read)]]


def test_qc_multiqc_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_dir = tmp_path / "fastqc"
    input_dir.mkdir()
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "multiqc" if name == "multiqc" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc(backend="multiqc", input_dir=input_dir, output_dir=tmp_path / "multiqc")

    assert calls == [["multiqc", str(input_dir), "--outdir", str(tmp_path / "multiqc")]]


def test_qc_multiqc_force_reaches_tool(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    input_dir = tmp_path / "fastqc"
    output_dir = tmp_path / "multiqc"
    input_dir.mkdir()
    output_dir.mkdir()
    (output_dir / "multiqc_report.html").write_text("old", encoding="utf-8")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "multiqc" if name == "multiqc" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc(backend="multiqc", input_dir=input_dir, output_dir=output_dir, force=True)

    assert calls == [["multiqc", str(input_dir), "--outdir", str(output_dir), "--force"]]


def test_qc_qiime2_demux_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc(backend="qiime2-demux", demux=demux, output=tmp_path / "demux.qzv")

    assert calls == [
        [
            "qiime",
            "demux",
            "summarize",
            "--i-data",
            str(demux),
            "--o-visualization",
            str(tmp_path / "demux.qzv"),
        ]
    ]


def test_qc_qiime2_demux_force_unlinks_existing_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    output = touch(tmp_path / "demux.qzv")

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        assert not output.exists()
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc(backend="qiime2-demux", demux=demux, output=output, force=True)


def test_trim_fastp_single_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "fastp" if name == "fastp" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    trim(
        backend="fastp",
        read1=read,
        output1=tmp_path / "trimmed_R1.fastq.gz",
        html=tmp_path / "fastp.html",
        json_report=tmp_path / "fastp.json",
        threads=4,
    )

    assert calls == [
        [
            "fastp",
            "--in1",
            str(read),
            "--out1",
            str(tmp_path / "trimmed_R1.fastq.gz"),
            "--html",
            str(tmp_path / "fastp.html"),
            "--json",
            str(tmp_path / "fastp.json"),
            "--thread",
            "4",
        ]
    ]


def test_trim_fastp_paired_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read1 = touch(tmp_path / "sample_R1.fastq.gz")
    read2 = touch(tmp_path / "sample_R2.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "fastp" if name == "fastp" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    trim(
        backend="fastp",
        read1=read1,
        read2=read2,
        output1=tmp_path / "trimmed_R1.fastq.gz",
        output2=tmp_path / "trimmed_R2.fastq.gz",
        html=tmp_path / "fastp.html",
        json_report=tmp_path / "fastp.json",
    )

    command = calls[0]
    assert command[:2] == ["fastp", "--in1"]
    assert "--in2" in command
    assert str(read2) in command
    assert "--out2" in command
    assert str(tmp_path / "trimmed_R2.fastq.gz") in command


def test_trim_planned_backend_message(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="registered but not implemented"):
        trim(backend="cutadapt", read1=read, output1=tmp_path / "trimmed.fastq.gz")


def test_cli_exposes_qc_trim_and_reports_missing_fastp(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    monkeypatch.setattr("shutil.which", lambda name: None)
    runner = CliRunner()

    methods = runner.invoke(app, ["methods"])
    assert methods.exit_code == 0
    assert "qc" in methods.stdout
    assert "trim" in methods.stdout

    result = runner.invoke(
        app,
        [
            "trim",
            "--backend",
            "fastp",
            "--read1",
            str(read),
            "--output1",
            str(tmp_path / "trimmed.fastq.gz"),
            "--html",
            str(tmp_path / "fastp.html"),
            "--json-report",
            str(tmp_path / "fastp.json"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None

    result = runner.invoke(
        app,
        [
            "qc",
            "--backend",
            "fastqc",
            "--input",
            str(read),
            "--output-dir",
            str(tmp_path / "qc"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
