from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata import validate_stage_result
from microsuite.methods.diversity_calc import diversity_calc
from microsuite.methods.tax_classify import tax_classify


def _envelope(run_dir: Path) -> dict:
    files = list((run_dir / "stage-results").glob("*.json"))
    assert len(files) == 1, files
    return json.loads(files[0].read_text())


def _fake_qiime(returncode: int):
    def fake_run(command, *, check, text, capture_output, **kw):
        return subprocess.CompletedProcess(command, returncode, "", "err" if returncode else "")

    return fake_run


def test_diversity_emits_stage_envelope(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr("subprocess.run", _fake_qiime(0))
    table = tmp_path / "table.qza"
    table.write_text("x")
    output = tmp_path / "shannon.qza"
    diversity_calc(backend="qiime2", metric="shannon", table=table, output=output)
    env = _envelope(tmp_path)
    assert env["stage"] == "diversity" and env["backend"] == "qiime2"
    assert env["status"] == "completed"
    assert env["subprocesses"] and env["subprocesses"][0]["status"] == "completed"
    assert env["outputs"][0]["kind"] == "alpha_diversity"
    assert env["metrics"] == {}  # qza metric extraction deferred
    assert validate_stage_result(env) == []


def test_taxonomy_emits_failed_envelope_on_nonzero(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    monkeypatch.setattr("subprocess.run", _fake_qiime(1))
    rep = tmp_path / "rep.qza"
    rep.write_text("x")
    clf = tmp_path / "classifier.qza"
    clf.write_text("c")
    output = tmp_path / "taxonomy.qza"
    with pytest.raises(MicrobiomeSuiteError):
        tax_classify(backend="qiime2", rep_seqs=rep, output=output, classifier=str(clf))
    env = _envelope(tmp_path)
    assert env["stage"] == "taxonomy" and env["backend"] == "qiime2"
    assert env["status"] == "failed"
    assert env["subprocesses"][0]["status"] == "failed"
    assert validate_stage_result(env) == []
