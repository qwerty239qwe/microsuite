from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

import microsuite.api as api
from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.decontam import decontam
from microsuite.methods.evaluate import evaluate
from microsuite.methods.qc_filter import qc_filter


def touch(path: Path) -> Path:
    path.write_text("placeholder", encoding="utf-8")
    return path


def qiime_only(name: str) -> str | None:
    return "qiime" if name == "qiime" else None


def test_python_sdk_facade_exports_qiime2_quality_control_methods() -> None:
    assert api.qc_filter is qc_filter
    assert api.decontam is decontam
    assert api.evaluate is evaluate
    assert "qc_filter" in api.__all__
    assert "decontam" in api.__all__
    assert "evaluate" in api.__all__


def test_qc_filter_qiime2_filter_reads_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    database = touch(tmp_path / "human-bowtie2-index.qza")
    output = tmp_path / "filtered.qza"
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", qiime_only)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc_filter(
        backend="qiime2-filter-reads",
        demux=demux,
        database=database,
        output=output,
        threads=4,
        mode="local",
        sensitivity="very-sensitive",
        exclude=True,
    )

    assert calls == [
        [
            "qiime",
            "quality-control",
            "filter-reads",
            "--i-demultiplexed-sequences",
            str(demux),
            "--i-database",
            str(database),
            "--p-n-threads",
            "4",
            "--p-mode",
            "local",
            "--p-sensitivity",
            "very-sensitive",
            "--p-exclude-seqs",
            "--o-filtered-sequences",
            str(output),
        ]
    ]


def test_qc_filter_qiime2_filter_reads_keep_matches_flag(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    database = touch(tmp_path / "human-bowtie2-index.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    qc_filter(
        backend="qiime2-filter-reads",
        demux=demux,
        database=database,
        output=tmp_path / "kept.qza",
        exclude=False,
    )

    assert "--p-no-exclude-seqs" in calls[0]
    assert "--p-exclude-seqs" not in calls[0]


def test_cli_qc_filter_qiime2_run_dir_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    database = touch(tmp_path / "human-bowtie2-index.qza")
    run_dir = tmp_path / "run"

    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    result = CliRunner().invoke(
        app,
        [
            "qc_filter",
            "--backend",
            "qiime2-filter-reads",
            "--demux",
            str(demux),
            "--database",
            str(database),
            "-o",
            str(tmp_path / "filtered.qza"),
            "--run-dir",
            str(run_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "qc_filter"
    assert run["backend"] == "qiime2-filter-reads"
    assert "quality-control" in run["command"]
    assert (run_dir / "command.txt").exists()
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    assert (run_dir / "stderr.log").read_text(encoding="utf-8") == ""
    assert "command_end" in (run_dir / "events.jsonl").read_text(encoding="utf-8")


def test_qc_filter_qiime2_exclude_seqs_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = touch(tmp_path / "rep-seqs.qza")
    reference = touch(tmp_path / "host-reference.qza")
    hits = tmp_path / "hits.qza"
    misses = tmp_path / "misses.qza"
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", qiime_only)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc_filter(
        backend="qiime2-exclude-seqs",
        query_sequences=query,
        reference_sequences=reference,
        sequence_hits=hits,
        sequence_misses=misses,
        method="vsearch",
        perc_identity=0.95,
        perc_query_aligned=0.9,
        threads=8,
    )

    assert calls == [
        [
            "qiime",
            "quality-control",
            "exclude-seqs",
            "--i-query-sequences",
            str(query),
            "--i-reference-sequences",
            str(reference),
            "--p-method",
            "vsearch",
            "--p-perc-identity",
            "0.95",
            "--p-perc-query-aligned",
            "0.9",
            "--p-threads",
            "8",
            "--o-sequence-hits",
            str(hits),
            "--o-sequence-misses",
            str(misses),
        ]
    ]


def test_qc_filter_qiime2_exclude_seqs_blast_omits_threads(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    query = touch(tmp_path / "rep-seqs.qza")
    reference = touch(tmp_path / "host-reference.qza")
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    qc_filter(
        backend="qiime2-exclude-seqs",
        query_sequences=query,
        reference_sequences=reference,
        sequence_hits=tmp_path / "hits.qza",
        sequence_misses=tmp_path / "misses.qza",
        method="blast",
        threads=8,
    )

    assert "--p-threads" not in calls[0]


def test_qc_filter_qiime2_bowtie2_build_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sequences = touch(tmp_path / "human-reference-seqs.qza")
    output = tmp_path / "human-bowtie2-index.qza"
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    qc_filter(
        backend="qiime2-bowtie2-build",
        sequences=sequences,
        output=output,
        threads=6,
    )

    assert calls == [
        [
            "qiime",
            "quality-control",
            "bowtie2-build",
            "--i-sequences",
            str(sequences),
            "--p-n-threads",
            "6",
            "--o-database",
            str(output),
        ]
    ]


def test_decontam_qiime2_identify_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = touch(tmp_path / "table.qza")
    metadata = touch(tmp_path / "metadata.tsv")
    output = tmp_path / "scores.qza"
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", qiime_only)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    decontam(
        backend="qiime2-decontam",
        table=table,
        metadata=metadata,
        output=output,
        method="prevalence",
        prev_control_column="sample_type",
        prev_control_indicator="blank",
    )

    assert calls == [
        [
            "qiime",
            "quality-control",
            "decontam-identify",
            "--i-table",
            str(table),
            "--m-metadata-file",
            str(metadata),
            "--p-method",
            "prevalence",
            "--p-prev-control-column",
            "sample_type",
            "--p-prev-control-indicator",
            "blank",
            "--o-decontam-scores",
            str(output),
        ]
    ]


def test_decontam_prevalence_documents_required_control_metadata(tmp_path: Path) -> None:
    table = touch(tmp_path / "table.qza")
    metadata = touch(tmp_path / "metadata.tsv")

    with pytest.raises(MicrobiomeSuiteError, match="prev-control-column"):
        decontam(
            backend="qiime2-decontam",
            table=table,
            metadata=metadata,
            output=tmp_path / "scores.qza",
        )


def test_qiime2_quality_control_missing_qiime_reports_task(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    database = touch(tmp_path / "database.qza")
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="filter-reads requires"):
        qc_filter(
            backend="qiime2-filter-reads",
            demux=demux,
            database=database,
            output=tmp_path / "filtered.qza",
        )


def test_evaluate_qiime2_taxonomy_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    expected = touch(tmp_path / "expected-taxonomy.qza")
    observed = touch(tmp_path / "observed-taxonomy.qza")
    table = touch(tmp_path / "table.qza")
    output = tmp_path / "taxonomy-evaluation.qzv"
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", qiime_only)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    evaluate(
        backend="qiime2-taxonomy",
        expected_taxa=expected,
        observed_taxa=observed,
        feature_table=table,
        output=output,
        depth=7,
    )

    assert calls == [
        [
            "qiime",
            "quality-control",
            "evaluate-taxonomy",
            "--i-expected-taxa",
            str(expected),
            "--i-observed-taxa",
            str(observed),
            "--i-feature-table",
            str(table),
            "--p-depth",
            "7",
            "--o-visualization",
            str(output),
        ]
    ]


def test_cli_exposes_qiime2_quality_control_method_groups(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    demux = touch(tmp_path / "demux.qza")
    database = touch(tmp_path / "database.qza")
    monkeypatch.setattr("shutil.which", lambda name: None)
    runner = CliRunner()

    methods = runner.invoke(app, ["methods"])
    assert methods.exit_code == 0
    assert "qc_filter" in methods.stdout
    assert "decontam" in methods.stdout
    assert "evaluate" in methods.stdout
    assert "qiime2-filter-reads" in methods.stdout
    assert "qiime2-bowtie2-build" in methods.stdout
    assert "qiime2-decontam" in methods.stdout
    assert "qiime2-taxonomy" in methods.stdout

    result = runner.invoke(
        app,
        [
            "qc_filter",
            "--backend",
            "qiime2-filter-reads",
            "--demux",
            str(demux),
            "--database",
            str(database),
            "--output",
            str(tmp_path / "filtered.qza"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
