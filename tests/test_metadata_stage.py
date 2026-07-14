from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata import (
    Artifact,
    ArtifactCount,
    ProvenanceFile,
    WorkflowContext,
    stage_execution,
    validate_stage_result,
)


def _envelopes(run_dir: Path) -> list[Path]:
    return sorted((run_dir / "stage-results").glob("*.json"))


def _one(run_dir: Path) -> tuple[dict, Path]:
    files = _envelopes(run_dir)
    assert len(files) == 1, files
    return json.loads(files[0].read_text()), files[0]


def test_success_writes_valid_enriched_envelope(tmp_path: Path) -> None:
    out = tmp_path / "table.tsv"
    out.write_text("x" * 10)
    with stage_execution(
        tmp_path,
        stage="denoise",
        backend="dada2-r",
        params={"trunc_len_f": 240, "auth_token": "s3cr3tvalue"},
        outputs=[
            Artifact(
                "table", out, format="tsv", kind="feature_table", count=ArtifactCount(5, "features")
            )
        ],
    ) as rec:
        rec.note_subprocess(
            ["Rscript", "denoise.R"], exit_code=0, duration_sec=1.0, status="completed"
        )
    env, path = _one(tmp_path)
    assert validate_stage_result(env) == []
    assert env["status"] == "completed" and env["exit_code"] == 0
    o = env["outputs"][0]
    assert (
        o["exists"] is True
        and o["bytes"] == 10
        and o["external"] is False
        and o["path"] == "table.tsv"
    )
    assert o["count"] == {"value": 5, "unit": "features"}
    assert env["params"]["auth_token"] == "***"
    assert env["metrics"] == {} and env["software"] == {} and env["reference_db"] is None
    assert path.name.startswith("denoise--dada2-r--attempt-1--stage-run-")
    assert not list((tmp_path / "stage-results").glob("*.tmp*"))


def test_directory_output_and_external_input(tmp_path: Path) -> None:
    d = tmp_path / "outdir"
    d.mkdir()
    ext_dir = tmp_path.parent / "ext_metadata"
    ext_dir.mkdir(exist_ok=True)
    ext = ext_dir / "in.txt"
    ext.write_text("y")
    with stage_execution(
        tmp_path,
        stage="denoise",
        backend="dada2-r",
        inputs=[Artifact("reads", ext, format="directory")],
        outputs=[Artifact("outdir", d, format="directory")],
    ):
        pass
    env, _ = _one(tmp_path)
    assert env["outputs"][0]["bytes"] is None and env["outputs"][0]["exists"] is True
    assert env["inputs"][0]["external"] is True
    assert Path(env["inputs"][0]["path"]).is_absolute()


def test_failure_records_and_reraises(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="boom"):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r") as rec:
            rec.note_subprocess(["r", "fail.R"], exit_code=1, duration_sec=1.0, status="failed")
            raise MicrobiomeSuiteError("boom")
    env, _ = _one(tmp_path)
    assert validate_stage_result(env) == []
    assert env["status"] == "failed" and env["exit_code"] == 1
    assert env["error"]["type"] == "MicrobiomeSuiteError"
    assert env["command"] == ["r", "fail.R"]


def test_timeout_status_null_exit(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r") as rec:
            rec.note_subprocess(["r"], exit_code=None, duration_sec=60.0, status="timed_out")
            raise MicrobiomeSuiteError("timed out")
    env, _ = _one(tmp_path)
    assert env["status"] == "timed_out" and env["exit_code"] is None
    assert validate_stage_result(env) == []


def test_launch_failed_preserves_command(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r") as rec:
            rec.note_subprocess(
                ["missing-bin"], exit_code=None, duration_sec=0.0, status="launch_failed"
            )
            raise FileNotFoundError("missing-bin")
    env, _ = _one(tmp_path)
    assert env["status"] == "failed"
    assert env["command"] == ["missing-bin"] and env["exit_code"] is None
    assert env["subprocesses"][0]["status"] == "launch_failed"
    assert validate_stage_result(env) == []


def test_multiple_subprocesses_one_envelope(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r") as rec:
            rec.note_subprocess(["a"], exit_code=0, duration_sec=1.0, status="completed")
            rec.note_subprocess(["b", "fail"], exit_code=2, duration_sec=1.0, status="failed")
            raise MicrobiomeSuiteError("x")
    env, _ = _one(tmp_path)
    assert len(env["subprocesses"]) == 2
    assert env["command"] == ["b", "fail"] and env["exit_code"] == 2


def test_retry_distinct_attempts_no_overwrite(tmp_path: Path) -> None:
    for _ in range(2):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r"):
            pass
    files = _envelopes(tmp_path)
    assert len(files) == 2
    attempts = {json.loads(f.read_text())["attempt"] for f in files}
    ids = {json.loads(f.read_text())["stage_run_id"] for f in files}
    assert attempts == {1, 2}
    assert len(ids) == 2


def test_task_defaults_to_stage(tmp_path: Path) -> None:
    with stage_execution(tmp_path, stage="denoise"):
        pass
    env, _ = _one(tmp_path)
    assert env["task"] == "denoise"


def test_success_path_invalid_payload_raises_and_diagnostic(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("microsuite.metadata.stage.validate_stage_result", lambda p: ["boom"])
    with pytest.raises(MicrobiomeSuiteError, match="invalid"):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r"):
            pass
    assert list((tmp_path / "stage-results" / "diagnostics").glob("*.invalid"))
    assert _envelopes(tmp_path) == []


def test_writer_failure_after_success_raises(tmp_path: Path, monkeypatch) -> None:
    def boom(*a, **k):
        raise OSError("disk full")

    monkeypatch.setattr("microsuite.metadata.stage._atomic_write", boom)
    with pytest.raises(MicrobiomeSuiteError, match="Failed to write"):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r"):
            pass


def test_explicit_workflow_context(tmp_path: Path) -> None:
    ctx = WorkflowContext(
        run_id="wf-run", workflow_id="wid", workflow_run_id="wrid", dataset_id="ds"
    )
    with stage_execution(tmp_path, stage="denoise", workflow_context=ctx):
        pass
    env, _ = _one(tmp_path)
    assert env["run_id"] == "wf-run"
    assert env["workflow_run_id"] == "wrid" and env["dataset_id"] == "ds"


def test_env_workflow_identity_and_standalone_fallback(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_WORKFLOW_RUN_ID", "env-wrid")
    with stage_execution(tmp_path, stage="denoise"):
        pass
    env, _ = _one(tmp_path)
    assert env["workflow_run_id"] == "env-wrid"
    assert env["run_id"].startswith("standalone-stage-run-")
    assert env["dataset_id"] is None


def test_cancellation(tmp_path: Path) -> None:
    with pytest.raises(KeyboardInterrupt):
        with stage_execution(tmp_path, stage="denoise", backend="dada2-r"):
            raise KeyboardInterrupt
    env, _ = _one(tmp_path)
    assert env["status"] == "cancelled" and env["exit_code"] is None
    assert env["error"]["type"] == "KeyboardInterrupt"
    assert validate_stage_result(env) == []


def test_run_dir_none_writes_nothing(tmp_path: Path) -> None:
    with stage_execution(None, stage="denoise") as rec:
        rec.note_subprocess(["r"], exit_code=0, duration_sec=1.0, status="completed")
    assert not (tmp_path / "stage-results").exists()


def test_provenance_reference_exists(tmp_path: Path) -> None:
    prov = tmp_path / "dada2_denoise_manifest.json"
    prov.write_text("{}")
    with stage_execution(
        tmp_path,
        stage="denoise",
        backend="dada2-r",
        provenance_files=[ProvenanceFile("dada2_manifest", prov)],
    ):
        pass
    env, _ = _one(tmp_path)
    assert env["provenance_files"][0]["exists"] is True
    assert env["provenance_files"][0]["path"] == "dada2_denoise_manifest.json"
