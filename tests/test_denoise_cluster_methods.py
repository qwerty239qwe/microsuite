from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

import microsuite.api as api
from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.cluster import cluster
from microsuite.methods.denoise import denoise


def touch(path: Path) -> Path:
    path.write_text("placeholder", encoding="utf-8")
    return path


def test_python_sdk_facade_exports_denoise() -> None:
    assert api.denoise is denoise
    assert "denoise" in api.__all__


def test_denoise_qiime2_dada2_single_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="qiime2-dada2",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        trunc_len=150,
        trim_left=5,
        threads=2,
    )

    assert calls == [
        [
            "qiime",
            "dada2",
            "denoise-single",
            "--i-demultiplexed-seqs",
            str(demux),
            "--p-trim-left",
            "5",
            "--p-trunc-len",
            "150",
            "--o-table",
            str(tmp_path / "table.qza"),
            "--o-representative-sequences",
            str(tmp_path / "rep-seqs.qza"),
            "--o-denoising-stats",
            str(tmp_path / "stats.qza"),
            "--p-n-threads",
            "2",
        ]
    ]


def test_denoise_qiime2_deblur_requires_positive_trim_length(tmp_path: Path) -> None:
    demux = touch(tmp_path / "demux.qza")

    with pytest.raises(MicrobiomeSuiteError, match="--trunc-len"):
        denoise(
            backend="qiime2-deblur",
            demux=demux,
            output_table=tmp_path / "table.qza",
            output_rep_seqs=tmp_path / "rep-seqs.qza",
            output_stats=tmp_path / "stats.qza",
        )


def test_cli_denoise_qiime2_run_dir_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    run_dir = tmp_path / "run"

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    result = CliRunner().invoke(
        app,
        [
            "denoise",
            "--backend",
            "qiime2-dada2",
            "--demux",
            str(demux),
            "--output-table",
            str(tmp_path / "table.qza"),
            "--output-rep-seqs",
            str(tmp_path / "rep-seqs.qza"),
            "--output-stats",
            str(tmp_path / "stats.qza"),
            "--trunc-len",
            "150",
            "--run-dir",
            str(run_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "denoise"
    assert run["backend"] == "qiime2-dada2"
    assert "dada2" in run["command"]
    assert (run_dir / "command.txt").exists()
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    assert (run_dir / "stderr.log").read_text(encoding="utf-8") == ""
    events = (run_dir / "events.jsonl").read_text(encoding="utf-8")
    assert "command_start" in events
    assert "command_end" in events


def test_denoise_qiime2_dada2_paired_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="qiime2-dada2",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        paired=True,
        trim_left_f=7,
        trunc_len_f=151,
        trim_left_r=11,
        trunc_len_r=149,
        threads=4,
    )

    command = calls[0]
    assert command[:3] == ["qiime", "dada2", "denoise-paired"]
    assert "--p-trim-left-f" in command
    assert "7" in command
    assert "--p-trunc-len-f" in command
    assert "151" in command
    assert "--p-trim-left-r" in command
    assert "11" in command
    assert "--p-trunc-len-r" in command
    assert "149" in command


def test_denoise_qiime2_deblur_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="qiime2-deblur",
        demux=demux,
        output_table=tmp_path / "table.qza",
        output_rep_seqs=tmp_path / "rep-seqs.qza",
        output_stats=tmp_path / "stats.qza",
        trim_left=3,
        trunc_len=120,
        threads=2,
    )

    assert calls == [
        [
            "qiime",
            "deblur",
            "denoise-16S",
            "--i-demultiplexed-seqs",
            str(demux),
            "--p-trim-length",
            "120",
            "--p-left-trim-len",
            "3",
            "--p-jobs-to-start",
            "2",
            "--o-table",
            str(tmp_path / "table.qza"),
            "--o-representative-sequences",
            str(tmp_path / "rep-seqs.qza"),
            "--o-stats",
            str(tmp_path / "stats.qza"),
        ]
    ]


def test_denoise_dada2_r_builds_rscript_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    denoise(
        backend="dada2-r",
        demux=reads,
        output_table=tmp_path / "table.tsv",
        output_rep_seqs=tmp_path / "rep-seqs.fasta",
        output_stats=tmp_path / "stats.tsv",
        paired=True,
        trim_left_f=7,
        trunc_len_f=151,
        trim_left_r=11,
        trunc_len_r=149,
        threads=4,
    )

    command = calls[0]
    assert command[0] == "Rscript"
    assert command[1].endswith(str(Path("microsuite/resources/dada2_denoise.R")))
    assert command[2] == "--input-dir"
    assert str(reads) in command
    assert "--paired" in command
    assert "--trunc-len-f" in command
    assert "151" in command
    assert "--threads" in command
    assert "4" in command


def test_denoise_dada2_r_missing_rscript(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="R/DADA2 denoising requires"):
        denoise(
            backend="dada2-r",
            demux=reads,
            output_table=tmp_path / "table.tsv",
            output_rep_seqs=tmp_path / "rep-seqs.fasta",
            output_stats=tmp_path / "stats.tsv",
        )


def test_cluster_vsearch_builds_qiime2_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = touch(tmp_path / "table.qza")
    rep_seqs = touch(tmp_path / "rep-seqs.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    cluster(
        backend="vsearch",
        table=table,
        rep_seqs=rep_seqs,
        output_table=tmp_path / "clustered-table.qza",
        output_rep_seqs=tmp_path / "clustered-rep-seqs.qza",
        identity=0.97,
    )

    assert calls == [
        [
            "qiime",
            "vsearch",
            "cluster-features-de-novo",
            "--i-table",
            str(table),
            "--i-sequences",
            str(rep_seqs),
            "--p-perc-identity",
            "0.97",
            "--o-clustered-table",
            str(tmp_path / "clustered-table.qza"),
            "--o-clustered-sequences",
            str(tmp_path / "clustered-rep-seqs.qza"),
        ]
    ]


def test_cluster_usearch_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    rep_seqs = touch(tmp_path / "rep-seqs.fasta")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "usearch" if name == "usearch" else None)

    def fake_run(
        command: list[str], *, check: bool, text: bool, capture_output: bool
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    cluster(
        backend="usearch",
        rep_seqs=rep_seqs,
        output_table=tmp_path / "clusters.uc",
        output_rep_seqs=tmp_path / "centroids.fasta",
        identity=0.99,
    )

    assert calls == [
        [
            "usearch",
            "-cluster_fast",
            str(rep_seqs),
            "-id",
            "0.99",
            "-centroids",
            str(tmp_path / "centroids.fasta"),
            "-uc",
            str(tmp_path / "clusters.uc"),
        ]
    ]


def test_cluster_usearch_missing_binary_reports_tool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = touch(tmp_path / "rep-seqs.fasta")
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="USEARCH clustering requires"):
        cluster(
            backend="usearch",
            rep_seqs=rep_seqs,
            output_table=tmp_path / "clusters.uc",
            output_rep_seqs=tmp_path / "centroids.fasta",
        )


def test_cli_exposes_denoise_cluster_and_reports_missing_qiime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    monkeypatch.setattr("shutil.which", lambda name: None)
    runner = CliRunner()

    methods = runner.invoke(app, ["methods"])
    assert methods.exit_code == 0
    assert "denoise" in methods.stdout
    assert "cluster" in methods.stdout
    assert "usearch" in methods.stdout

    result = runner.invoke(
        app,
        [
            "denoise",
            "--backend",
            "qiime2-dada2",
            "--demux",
            str(demux),
            "--output-table",
            str(tmp_path / "table.qza"),
            "--output-rep-seqs",
            str(tmp_path / "rep-seqs.qza"),
            "--output-stats",
            str(tmp_path / "stats.qza"),
            "--trunc-len",
            "150",
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None

    table = touch(tmp_path / "table.qza")
    rep_seqs = touch(tmp_path / "rep-seqs.qza")
    result = runner.invoke(
        app,
        [
            "cluster",
            "--backend",
            "vsearch",
            "--table",
            str(table),
            "--rep-seqs",
            str(rep_seqs),
            "--output-table",
            str(tmp_path / "clustered-table.qza"),
            "--output-rep-seqs",
            str(tmp_path / "clustered-rep-seqs.qza"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
