from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata import stage_execution, validate_stage_result
from microsuite.runtime.runner import run_command


def _envelope(run_dir: Path) -> dict:
    files = list((run_dir / "stage-results").glob("*.json"))
    assert len(files) == 1, files
    return json.loads(files[0].read_text())


def test_run_command_contributes_completed_subprocess(tmp_path: Path) -> None:
    with stage_execution(tmp_path, stage="probe", backend="py"):
        run_command([sys.executable, "-c", "pass"], "failed", run_dir=tmp_path / "cmd")
    env = _envelope(tmp_path)
    assert env["status"] == "completed"
    assert len(env["subprocesses"]) == 1
    assert env["subprocesses"][0]["status"] == "completed"
    assert env["subprocesses"][0]["exit_code"] == 0


def test_run_command_contributes_independent_of_run_dir(tmp_path: Path) -> None:
    with stage_execution(tmp_path, stage="probe", backend="py"):
        run_command([sys.executable, "-c", "pass"], "failed", run_dir=None)
    env = _envelope(tmp_path)
    assert len(env["subprocesses"]) == 1 and env["subprocesses"][0]["status"] == "completed"


def test_run_command_nonzero_marks_failed_subprocess_and_stage(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        with stage_execution(tmp_path, stage="probe", backend="py"):
            run_command([sys.executable, "-c", "import sys; sys.exit(3)"], "boom", run_dir=None)
    env = _envelope(tmp_path)
    assert env["status"] == "failed"
    assert env["subprocesses"][0]["status"] == "failed"
    assert env["command"] == [sys.executable, "-c", "import sys; sys.exit(3)"]
    assert env["exit_code"] == 3
    assert validate_stage_result(env) == []


def test_run_command_launch_failure_records_launch_failed(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        with stage_execution(tmp_path, stage="probe", backend="none"):
            run_command(["definitely-not-a-real-binary-xyz"], "boom", run_dir=None)
    env = _envelope(tmp_path)
    assert env["subprocesses"][0]["status"] == "launch_failed"
    assert env["subprocesses"][0]["exit_code"] is None
    assert env["status"] == "failed"
    assert validate_stage_result(env) == []


def test_no_active_stage_writes_no_envelope_but_keeps_results_manifest(tmp_path: Path) -> None:
    run_command([sys.executable, "-c", "pass"], "failed", run_dir=tmp_path)
    assert not (tmp_path / "stage-results").exists()
    assert (tmp_path / "microsuite-results.json").exists()
