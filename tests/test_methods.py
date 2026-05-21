from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.diversity_calc import diversity_calc
from microsuite.methods.tax_classify import tax_classify


def test_tax_classify_qiime2_requires_classifier(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.qza"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--classifier"):
        tax_classify(backend="qiime2", rep_seqs=rep_seqs, output=tmp_path / "taxonomy.qza")


def test_tax_classify_planned_method_message(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.fastq"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="registered but not implemented"):
        tax_classify(backend="kraken2", rep_seqs=rep_seqs, output=tmp_path / "taxonomy.tsv")


def test_diversity_calc_qiime2_requires_phylogeny_for_unifrac(tmp_path: Path) -> None:
    table = tmp_path / "table.qza"
    table.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--phylogeny"):
        diversity_calc(
            backend="qiime2",
            metric="weighted-unifrac",
            table=table,
            output=tmp_path / "weighted-unifrac.qza",
        )


def test_diversity_calc_qiime2_reports_missing_qiime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.qza"
    table.write_text("placeholder", encoding="utf-8")

    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="qiime"):
        diversity_calc(
            backend="qiime2",
            metric="bray-curtis",
            table=table,
            output=tmp_path / "bray-curtis.qza",
        )


def test_cli_tax_classify_qiime2_run_dir_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = tmp_path / "rep-seqs.qza"
    classifier = tmp_path / "classifier.qza"
    rep_seqs.write_text("placeholder", encoding="utf-8")
    classifier.write_text("placeholder", encoding="utf-8")
    run_dir = tmp_path / "run"

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    result = CliRunner().invoke(
        app,
        [
            "tax_classify",
            "--backend",
            "qiime2",
            "--rep-seqs",
            str(rep_seqs),
            "--classifier",
            str(classifier),
            "-o",
            str(tmp_path / "taxonomy.qza"),
            "--run-dir",
            str(run_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "tax_classify"
    assert run["backend"] == "qiime2"
    assert "feature-classifier" in run["command"]
    assert (run_dir / "command.txt").exists()
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    assert (run_dir / "stderr.log").read_text(encoding="utf-8") == ""
    assert "command_end" in (run_dir / "events.jsonl").read_text(encoding="utf-8")


def test_cli_diversity_calc_qiime2_run_dir_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.qza"
    table.write_text("placeholder", encoding="utf-8")
    run_dir = tmp_path / "run"

    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    result = CliRunner().invoke(
        app,
        [
            "diversity_calc",
            "--backend",
            "qiime2",
            "--metric",
            "bray-curtis",
            "--table",
            str(table),
            "-o",
            str(tmp_path / "bray-curtis.qza"),
            "--run-dir",
            str(run_dir),
        ],
    )

    assert result.exit_code == 0, result.stdout
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "diversity_calc"
    assert run["backend"] == "qiime2"
    assert "diversity-lib" in run["command"]
    assert (run_dir / "command.txt").exists()
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    assert (run_dir / "stderr.log").read_text(encoding="utf-8") == ""
    assert "command_end" in (run_dir / "events.jsonl").read_text(encoding="utf-8")
