from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.qiime2_wrappers import (
    demux,
    feature_summarize,
    metadata_tabulate,
    phylogeny,
    tax_train,
)
from microsuite.workflows.moving_pictures_qiime2 import run_moving_pictures_qiime2


def touch(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("placeholder", encoding="utf-8")
    return path


def qiime_only(name: str) -> str | None:
    return "qiime" if name == "qiime" else None


def fake_success(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
    return subprocess.CompletedProcess(command, 0, "ok\n", "")


def test_metadata_tabulate_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = touch(tmp_path / "metadata.tsv")
    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: calls.append(command) or fake_success(command),
    )

    metadata_tabulate(
        backend="qiime2",
        input_file=metadata,
        output=tmp_path / "metadata.qzv",
        run_dir=tmp_path / "run",
    )

    assert calls == [
        [
            "qiime",
            "metadata",
            "tabulate",
            "--m-input-file",
            str(metadata),
            "--o-visualization",
            str(tmp_path / "metadata.qzv"),
        ]
    ]
    run = json.loads((tmp_path / "run" / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "metadata_tabulate"
    assert run["backend"] == "qiime2"


def test_new_wrapper_errors_are_clear(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="Unsupported metadata_tabulate backend"):
        metadata_tabulate(backend="native", input_file=None, output=tmp_path / "x.qzv")
    with pytest.raises(MicrobiomeSuiteError, match="--input-file is required"):
        metadata_tabulate(backend="qiime2", input_file=None, output=tmp_path / "x.qzv")
    with pytest.raises(MicrobiomeSuiteError, match="metadata tabulate requires"):
        metadata_tabulate(
            backend="qiime2",
            input_file=touch(tmp_path / "metadata.tsv"),
            output=tmp_path / "x.qzv",
        )


def test_demux_and_feature_summarize_commands(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    seqs = touch(tmp_path / "emp.qza")
    metadata = touch(tmp_path / "metadata.tsv")
    table = touch(tmp_path / "table.qza")
    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: calls.append(command) or fake_success(command),
    )

    demux(
        backend="qiime2-emp-single",
        seqs=seqs,
        metadata=metadata,
        barcode_column="barcode-sequence",
        output_demux=tmp_path / "demux.qza",
        output_details=tmp_path / "details.qza",
    )
    feature_summarize(
        backend="qiime2",
        mode="summarize",
        table=table,
        metadata=metadata,
        output=tmp_path / "table.qzv",
    )

    assert calls[0][1:3] == ["demux", "emp-single"]
    assert "--m-barcodes-column" in calls[0]
    assert calls[1][1:3] == ["feature-table", "summarize"]
    assert "--m-sample-metadata-file" in calls[1]


def test_phylogeny_and_tax_train_commands(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    rep_seqs = touch(tmp_path / "rep-seqs.qza")
    ref_seqs = touch(tmp_path / "85_otus.qza")
    ref_taxonomy = touch(tmp_path / "ref-taxonomy.qza")
    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: calls.append(command) or fake_success(command),
    )

    phylogeny(
        backend="qiime2-mafft-fasttree",
        rep_seqs=rep_seqs,
        output_aligned=tmp_path / "aligned.qza",
        output_masked=tmp_path / "masked.qza",
        output_tree=tmp_path / "tree.qza",
        output_rooted_tree=tmp_path / "rooted-tree.qza",
        threads=2,
    )
    tax_train(
        backend="qiime2-naive-bayes",
        ref_seqs=ref_seqs,
        ref_taxonomy=ref_taxonomy,
        f_primer="GTGCCAGCMGCCGCGGTAA",
        r_primer="GGACTACHVGGGTWTCTAAT",
        trunc_len=120,
        output=tmp_path / "classifier.qza",
        run_dir=tmp_path / "tax-train-run",
    )

    assert calls[0][1:3] == ["phylogeny", "align-to-tree-mafft-fasttree"]
    assert calls[1][1:3] == ["feature-classifier", "extract-reads"]
    assert calls[2][1:3] == ["feature-classifier", "fit-classifier-naive-bayes"]
    assert (tmp_path / "tax-train-run" / "extract-reads" / "run.json").exists()
    assert (tmp_path / "tax-train-run" / "fit-classifier-naive-bayes" / "run.json").exists()


def test_cli_new_method_help_and_runtime(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    metadata = touch(tmp_path / "metadata.tsv")
    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr("subprocess.run", fake_success)
    runner = CliRunner()

    help_result = runner.invoke(app, ["metadata_tabulate", "--help"], terminal_width=160)
    assert help_result.exit_code == 0

    result = runner.invoke(
        app,
        [
            "metadata_tabulate",
            "--backend",
            "qiime2",
            "--input-file",
            str(metadata),
            "--output",
            str(tmp_path / "metadata.qzv"),
            "--run-dir",
            str(tmp_path / "cli-run"),
        ],
    )
    assert result.exit_code == 0, result.stdout
    run = json.loads((tmp_path / "cli-run" / "run.json").read_text(encoding="utf-8"))
    assert run["command"][1:3] == ["metadata", "tabulate"]


def test_workflow_moving_pictures_qiime2_static(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fake_fetch(output: Path, *, force: bool = False) -> None:
        output.mkdir(parents=True, exist_ok=True)
        touch(output / "sample-metadata.tsv")
        (output / "emp-single-end-sequences").mkdir(exist_ok=True)
        touch(output / "emp-single-end-sequences" / "sequences.fastq.gz")
        touch(output / "85_otus.qza")
        touch(output / "ref-taxonomy.qza")

    calls: list[list[str]] = []
    monkeypatch.setattr("shutil.which", qiime_only)
    monkeypatch.setattr(
        "microsuite.workflows.moving_pictures_qiime2.fetch_moving_pictures",
        fake_fetch,
    )
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: calls.append(command) or fake_success(command),
    )

    run_moving_pictures_qiime2(output=tmp_path / "run", force=True, threads=1)

    run = json.loads((tmp_path / "run" / "run.json").read_text(encoding="utf-8"))
    step_names = [step["name"] for step in run["steps"]]
    assert step_names[:5] == [
        "metadata_tabulate",
        "import_emp_single_end",
        "demux_emp_single",
        "demux_summarize",
        "dada2_denoise_single",
    ]
    assert "fit_classifier_naive_bayes" in step_names
    assert "ancombc_subject_level_6" in step_names
    assert len(calls) == len(step_names)
    assert (tmp_path / "run" / "runtime" / "01-metadata_tabulate" / "run.json").exists()
    assert (tmp_path / "run" / "report.html").exists()


def test_workflow_reports_missing_qiime(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="requires the external 'qiime'"):
        run_moving_pictures_qiime2(output=tmp_path / "run")
