from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

import microsuite.api as api
from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.qc import qc
from microsuite.methods.trim import trim


def touch(path: Path) -> Path:
    path.write_text("placeholder", encoding="utf-8")
    return path


def test_python_sdk_facade_exports_trim() -> None:
    assert api.trim is trim
    assert "trim" in api.__all__


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


def test_qc_fastqc_auto_threads_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("os.cpu_count", lambda: 8)
    monkeypatch.setattr("shutil.which", lambda name: "fastqc" if name == "fastqc" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    qc(backend="fastqc", inputs=[read], output_dir=tmp_path / "qc", threads="auto")

    assert calls == [["fastqc", "--outdir", str(tmp_path / "qc"), "--threads", "7", str(read)]]


def test_cli_qc_fastqc_run_dir_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    run_dir = tmp_path / "run"

    monkeypatch.setattr("os.cpu_count", lambda: 8)
    monkeypatch.setattr("shutil.which", lambda name: "fastqc" if name == "fastqc" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    result = CliRunner().invoke(
        app,
        [
            "qc",
            "--backend",
            "fastqc",
            "--input",
            str(read),
            "--output-dir",
            str(tmp_path / "qc"),
            "--threads",
            "auto",
            "--run-dir",
            str(run_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    command_text = (run_dir / "command.txt").read_text(encoding="utf-8")
    assert command_text.startswith("fastqc --outdir ")
    assert "--threads 7" in command_text
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "qc"
    assert run["backend"] == "fastqc"


def test_cli_trim_cutadapt_run_dir_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    run_dir = tmp_path / "run"

    monkeypatch.setattr("shutil.which", lambda name: "cutadapt" if name == "cutadapt" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    result = CliRunner().invoke(
        app,
        [
            "trim",
            "--backend",
            "cutadapt",
            "--read1",
            str(read),
            "--output1",
            str(tmp_path / "trimmed_R1.fastq.gz"),
            "--adapter",
            "AGATCGGAAGAGC",
            "--run-dir",
            str(run_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "trim"
    assert run["backend"] == "cutadapt"
    assert "cutadapt" in run["command"]
    assert (run_dir / "command.txt").exists()
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    assert (run_dir / "stderr.log").read_text(encoding="utf-8") == ""
    events = (run_dir / "events.jsonl").read_text(encoding="utf-8")
    assert "command_start" in events
    assert "command_end" in events


def test_qc_fastqc_extract_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "fastqc" if name == "fastqc" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc(backend="fastqc", inputs=[read], output_dir=tmp_path / "qc", threads=2, extract=True)

    assert calls == [
        ["fastqc", "--outdir", str(tmp_path / "qc"), "--threads", "2", "--extract", str(read)]
    ]


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


def test_trim_cutadapt_single_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "cutadapt" if name == "cutadapt" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    trim(
        backend="cutadapt",
        read1=read,
        output1=tmp_path / "trimmed_R1.fastq.gz",
        adapter="AGATCGGAAGAGC",
        quality_cutoff="20",
        minimum_length="100",
        json_report=tmp_path / "cutadapt.json",
        threads=4,
    )

    assert calls == [
        [
            "cutadapt",
            "-a",
            "AGATCGGAAGAGC",
            "-q",
            "20",
            "-m",
            "100",
            "--json",
            str(tmp_path / "cutadapt.json"),
            "-j",
            "4",
            "-o",
            str(tmp_path / "trimmed_R1.fastq.gz"),
            str(read),
        ]
    ]


def test_trim_cutadapt_paired_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read1 = touch(tmp_path / "sample_R1.fastq.gz")
    read2 = touch(tmp_path / "sample_R2.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "cutadapt" if name == "cutadapt" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    trim(
        backend="cutadapt",
        read1=read1,
        read2=read2,
        output1=tmp_path / "trimmed_R1.fastq.gz",
        output2=tmp_path / "trimmed_R2.fastq.gz",
        front="^GTGYCAGCMGCCGCGGTAA",
        adapter2="AGATCGGAAGAGC",
        max_n="0",
        discard_untrimmed=True,
    )

    assert calls == [
        [
            "cutadapt",
            "-g",
            "^GTGYCAGCMGCCGCGGTAA",
            "-A",
            "AGATCGGAAGAGC",
            "--max-n",
            "0",
            "--discard-untrimmed",
            "-j",
            "1",
            "-o",
            str(tmp_path / "trimmed_R1.fastq.gz"),
            "-p",
            str(tmp_path / "trimmed_R2.fastq.gz"),
            str(read1),
            str(read2),
        ]
    ]


def test_trim_trimmomatic_single_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "shutil.which", lambda name: "trimmomatic" if name == "trimmomatic" else None
    )

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    trim(
        backend="trimmomatic",
        read1=read,
        output1=tmp_path / "trimmed_R1.fastq.gz",
        trimmomatic_steps=[
            "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10",
            "SLIDINGWINDOW:4:20",
            "MINLEN:100",
        ],
        threads=4,
    )

    assert calls == [
        [
            "trimmomatic",
            "SE",
            "-threads",
            "4",
            str(read),
            str(tmp_path / "trimmed_R1.fastq.gz"),
            "ILLUMINACLIP:TruSeq3-SE.fa:2:30:10",
            "SLIDINGWINDOW:4:20",
            "MINLEN:100",
        ]
    ]


def test_trim_trimmomatic_paired_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read1 = touch(tmp_path / "sample_R1.fastq.gz")
    read2 = touch(tmp_path / "sample_R2.fastq.gz")
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "shutil.which", lambda name: "trimmomatic" if name == "trimmomatic" else None
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, *, check, text, capture_output: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    trim(
        backend="trimmomatic",
        read1=read1,
        read2=read2,
        output1=tmp_path / "paired_R1.fastq.gz",
        output2=tmp_path / "paired_R2.fastq.gz",
        unpaired1=tmp_path / "unpaired_R1.fastq.gz",
        unpaired2=tmp_path / "unpaired_R2.fastq.gz",
        trimmomatic_steps=["LEADING:3", "TRAILING:3", "MINLEN:100"],
    )

    assert calls == [
        [
            "trimmomatic",
            "PE",
            "-threads",
            "1",
            str(read1),
            str(read2),
            str(tmp_path / "paired_R1.fastq.gz"),
            str(tmp_path / "unpaired_R1.fastq.gz"),
            str(tmp_path / "paired_R2.fastq.gz"),
            str(tmp_path / "unpaired_R2.fastq.gz"),
            "LEADING:3",
            "TRAILING:3",
            "MINLEN:100",
        ]
    ]


def test_trim_trim_galore_single_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "shutil.which", lambda name: "trim_galore" if name == "trim_galore" else None
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, *, check, text, capture_output: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    trim(
        backend="trim-galore",
        read1=read,
        output1=tmp_path / "trim_galore" / "sample_trimmed.fq.gz",
        quality_cutoff="20",
        minimum_length="100",
        basename="sample_trimmed",
        threads=2,
    )

    assert calls == [
        [
            "trim_galore",
            "--quality",
            "20",
            "--length",
            "100",
            "--cores",
            "2",
            "--output_dir",
            str(tmp_path / "trim_galore"),
            "--basename",
            "sample_trimmed",
            str(read),
        ]
    ]


def test_trim_trim_galore_v2_version_flag_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "shutil.which", lambda name: "trim_galore" if name == "trim_galore" else None
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, *, check, text, capture_output: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    trim(
        backend="trim-galore",
        read1=read,
        output1=tmp_path / "trim_galore" / "sample_trimmed.fq.gz",
        basename="sample_trimmed",
        trim_galore_version="v2",
        threads=4,
    )

    assert calls == [
        [
            "trim_galore",
            "--cores",
            "4",
            "--output_dir",
            str(tmp_path / "trim_galore"),
            "--basename",
            "sample_trimmed",
            "--engine",
            "v2",
            str(read),
        ]
    ]


def test_trim_trim_galore_rejects_unknown_version(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="trim-galore-version"):
        trim(
            backend="trim-galore",
            read1=read,
            output1=tmp_path / "sample_trimmed.fq.gz",
            trim_galore_version="classic",
        )


def test_trim_trim_galore_requires_output_basename_match(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="expected output1"):
        trim(
            backend="trim-galore",
            read1=read,
            output1=tmp_path / "custom.fastq.gz",
            basename="sample_trimmed",
        )


def test_trim_trimmomatic_rejects_cutadapt_options(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="not supported by --backend trimmomatic"):
        trim(
            backend="trimmomatic",
            read1=read,
            output1=tmp_path / "trimmed.fastq.gz",
            adapter="AGATCGGAAGAGC",
            trimmomatic_steps=["MINLEN:100"],
        )


def test_trim_galore_rejects_unsupported_shared_options(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="not supported by --backend trim-galore"):
        trim(
            backend="trim-galore",
            read1=read,
            output1=tmp_path / "sample_trimmed.fq.gz",
            front="^GTGYCAGCMGCCGCGGTAA",
        )


def test_trim_trim_galore_paired_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read1 = touch(tmp_path / "sample_R1.fastq.gz")
    read2 = touch(tmp_path / "sample_R2.fastq.gz")
    calls: list[list[str]] = []
    monkeypatch.setattr(
        "shutil.which", lambda name: "trim_galore" if name == "trim_galore" else None
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, *, check, text, capture_output: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    trim(
        backend="trim-galore",
        read1=read1,
        read2=read2,
        output1=tmp_path / "trim_galore" / "sample_R1_val_1.fq.gz",
        output2=tmp_path / "trim_galore" / "sample_R2_val_2.fq.gz",
        adapter="AGATCGGAAGAGC",
        adapter2="AGATCGGAAGAGC",
        threads=2,
    )

    assert calls == [
        [
            "trim_galore",
            "--paired",
            "--adapter",
            "AGATCGGAAGAGC",
            "--adapter2",
            "AGATCGGAAGAGC",
            "--cores",
            "2",
            "--output_dir",
            str(tmp_path / "trim_galore"),
            str(read1),
            str(read2),
        ]
    ]


def test_trim_cutadapt_requires_adapter_or_filter(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="adapter or filtering option"):
        trim(backend="cutadapt", read1=read, output1=tmp_path / "trimmed.fastq.gz")


def test_trim_cutadapt_discard_untrimmed_requires_adapter(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="requires at least one adapter"):
        trim(
            backend="cutadapt",
            read1=read,
            output1=tmp_path / "trimmed.fastq.gz",
            discard_untrimmed=True,
        )


def test_trim_cutadapt_missing_binary_reports_tool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="Cutadapt trimming requires"):
        trim(
            backend="cutadapt",
            read1=read,
            output1=tmp_path / "trimmed.fastq.gz",
            adapter="AGATCGGAAGAGC",
        )


def test_trim_cutadapt_rejects_html_report(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="--html is only supported"):
        trim(
            backend="cutadapt",
            read1=read,
            output1=tmp_path / "trimmed.fastq.gz",
            adapter="AGATCGGAAGAGC",
            html=tmp_path / "cutadapt.html",
        )


def test_trim_fastp_rejects_cutadapt_options(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="Cutadapt-specific trim options"):
        trim(
            backend="fastp",
            read1=read,
            output1=tmp_path / "trimmed.fastq.gz",
            adapter="AGATCGGAAGAGC",
        )


def test_trim_cutadapt_r2_adapter_requires_paired_input(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="R2 adapter options require"):
        trim(
            backend="cutadapt",
            read1=read,
            output1=tmp_path / "trimmed.fastq.gz",
            adapter2="AGATCGGAAGAGC",
        )


def test_trim_qiime2_cutadapt_planned_backend_message(tmp_path: Path) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")

    with pytest.raises(MicrobiomeSuiteError, match="registered but not implemented"):
        trim(backend="qiime2-cutadapt", read1=read, output1=tmp_path / "trimmed.fastq.gz")


def test_cli_trim_cutadapt_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read = touch(tmp_path / "sample_R1.fastq.gz")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "cutadapt" if name == "cutadapt" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    result = CliRunner().invoke(
        app,
        [
            "trim",
            "--backend",
            "cutadapt",
            "--read1",
            str(read),
            "--output1",
            str(tmp_path / "trimmed.fastq.gz"),
            "--adapter",
            "AGATCGGAAGAGC",
            "--quality-cutoff",
            "20",
            "--minimum-length",
            "100",
            "--json-report",
            str(tmp_path / "cutadapt.json"),
            "--threads",
            "4",
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert calls == [
        [
            "cutadapt",
            "-a",
            "AGATCGGAAGAGC",
            "-q",
            "20",
            "-m",
            "100",
            "--json",
            str(tmp_path / "cutadapt.json"),
            "-j",
            "4",
            "-o",
            str(tmp_path / "trimmed.fastq.gz"),
            str(read),
        ]
    ]


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
