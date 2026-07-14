from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata import validate_stage_result
from microsuite.methods.denoise import denoise


def _envelope(run_dir: Path) -> dict:
    files = list((run_dir / "stage-results").glob("*.json"))
    assert len(files) == 1, files
    return json.loads(files[0].read_text())


def _run_denoise(tmp_path: Path, monkeypatch, returncode: int) -> None:
    reads = tmp_path / "reads"
    reads.mkdir()
    out_table = tmp_path / "table.tsv"
    out_rep = tmp_path / "rep-seqs.fasta"
    out_stats = tmp_path / "stats.tsv"
    monkeypatch.setattr("shutil.which", lambda name: "Rscript" if name == "Rscript" else None)

    def fake_run(command, *, check, text, capture_output):
        if returncode == 0:
            out_table.write_text("#OTU\ts1\nASV1\t5\n")
            out_rep.write_text(">ASV1\nACGT\n")
            out_stats.write_text("\tinput\tfiltered\tmerged\tnonchim\ns1\t100\t90\t80\t70\n")
            (out_stats.parent / "dada2_r_params.json").write_text(
                json.dumps({"dada2_version": "1.30.0", "r_version": "4.3.1"})
            )
        return subprocess.CompletedProcess(command, returncode, "", "R error" if returncode else "")

    monkeypatch.setattr("subprocess.run", fake_run)
    denoise(
        backend="dada2-r",
        demux=reads,
        output_table=out_table,
        output_rep_seqs=out_rep,
        output_stats=out_stats,
        validate=False,
    )


def test_denoise_stage_envelope_finalized_after_provenance(tmp_path: Path, monkeypatch) -> None:
    _run_denoise(tmp_path, monkeypatch, returncode=0)
    env = _envelope(tmp_path)
    assert validate_stage_result(env) == []
    assert env["stage"] == "denoise" and env["backend"] == "dada2-r"
    assert env["status"] == "completed"
    assert env["subprocesses"][0]["status"] == "completed"
    outs = {o["kind"]: o for o in env["outputs"]}
    assert outs["feature_table"]["exists"] is True and outs["feature_table"]["path"] == "table.tsv"
    prov = {p["kind"]: p for p in env["provenance_files"]}
    assert prov["dada2_manifest"]["exists"] is True
    # the manifest really was written before finalization
    assert (tmp_path / "dada2_denoise_manifest.json").exists()
    # software versions lifted from the R params file (B)
    assert env["software"]["dada2"] == {"version": "1.30.0"}
    assert env["software"]["R"] == {"version": "4.3.1"}
    assert env["software"]["microsuite"]["version"]
    # retention metrics with units (C)
    m = env["metrics"]
    assert m["input_reads"] == {"value": 100, "unit": "reads"}
    assert m["nonchimeric_fraction"] == {"value": 0.7, "unit": "fraction"}
    assert m["nonchimeric_reads"] == {"value": 70, "unit": "reads"}


def test_denoise_stage_python_or_subprocess_failure(tmp_path: Path, monkeypatch) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        _run_denoise(tmp_path, monkeypatch, returncode=1)
    env = _envelope(tmp_path)
    assert env["status"] == "failed"
    assert env["subprocesses"][0]["status"] == "failed"
    assert env["exit_code"] == 1
    assert env["error"]["type"] == "MicrobiomeSuiteError"
    assert validate_stage_result(env) == []
