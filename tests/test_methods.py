from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.diversity_calc import diversity_calc
from microsuite.methods.qiime2_wrappers import phylogeny
from microsuite.methods.tax_classify import tax_classify


def test_tax_classify_qiime2_requires_classifier(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.qza"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--classifier"):
        tax_classify(backend="qiime2", rep_seqs=rep_seqs, output=tmp_path / "taxonomy.qza")


def test_tax_classify_kraken2_builds_command(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = tmp_path / "rep-seqs.fastq"
    database = tmp_path / "kraken-db"
    rep_seqs.write_text("placeholder", encoding="utf-8")
    database.mkdir()
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "kraken2" if name == "kraken2" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "", "")
        ),
    )

    tax_classify(
        backend="kraken2",
        rep_seqs=rep_seqs,
        classifier=database,
        output=tmp_path / "kraken-report.tsv",
        threads=3,
    )

    assert calls == [
        [
            "kraken2",
            "--db",
            str(database),
            "--threads",
            "3",
            "--report",
            str(tmp_path / "kraken-report.tsv"),
            "--output",
            str(tmp_path / "kraken-report.kraken"),
            str(rep_seqs),
        ]
    ]


def test_tax_classify_kraken2_requires_database(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.fastq"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--classifier"):
        tax_classify(backend="kraken2", rep_seqs=rep_seqs, output=tmp_path / "taxonomy.tsv")


def test_phylogeny_mafft_fasttree_builds_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = tmp_path / "rep-seqs.fasta"
    rep_seqs.write_text(">s1\nACGT\n>s2\nACGA\n", encoding="utf-8")
    calls: list[list[str]] = []

    def fake_which(name: str) -> str | None:
        return name if name in {"mafft", "FastTree"} else None

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        if command[0] == "mafft":
            return subprocess.CompletedProcess(command, 0, ">s1\nACGT\n>s2\nACGA\n", "")
        return subprocess.CompletedProcess(command, 0, "(s1:0.1,s2:0.1);\n", "")

    monkeypatch.setattr("shutil.which", fake_which)
    monkeypatch.setattr("subprocess.run", fake_run)

    phylogeny(
        backend="mafft-fasttree",
        rep_seqs=rep_seqs,
        output_aligned=tmp_path / "aligned.fasta",
        output_masked=tmp_path / "masked.fasta",
        output_tree=tmp_path / "tree.nwk",
        output_rooted_tree=tmp_path / "rooted-tree.nwk",
        threads=2,
        run_dir=tmp_path / "run",
    )

    assert calls[0] == ["mafft", "--auto", "--thread", "2", str(rep_seqs)]
    assert calls[1] == ["FastTree", "-nt", str(tmp_path / "aligned.fasta")]
    assert (tmp_path / "aligned.fasta").read_text(encoding="utf-8").startswith(">s1")
    assert (tmp_path / "tree.nwk").read_text(encoding="utf-8") == "(s1:0.1,s2:0.1);\n"
    assert (tmp_path / "run" / "mafft" / "run.json").exists()
    assert (tmp_path / "run" / "fasttree" / "run.json").exists()


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
