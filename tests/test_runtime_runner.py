from __future__ import annotations

import json
import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command


def test_run_command_writes_structured_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[list[str]] = []

    def fake_run(
        command: list[str],
        *,
        check: bool,
        text: bool,
        capture_output: bool,
    ) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "hello\n", "warn\n")

    monkeypatch.setattr("subprocess.run", fake_run)

    result = run_command(
        ["tool", "--flag"],
        "tool failed",
        run_dir=tmp_path,
        log=CommandLog(task="trim", backend="tool"),
    )

    assert result.returncode == 0
    assert calls == [["tool", "--flag"]]
    assert (tmp_path / "command.txt").read_text(encoding="utf-8") == "tool --flag\n"
    assert (tmp_path / "stdout.log").read_text(encoding="utf-8") == "hello\n"
    assert (tmp_path / "stderr.log").read_text(encoding="utf-8") == "warn\n"
    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert [event["event"] for event in events] == ["command_start", "command_end"]
    assert events[0]["task"] == "trim"
    assert events[0]["backend"] == "tool"
    assert events[1]["exit_code"] == 0
    run = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "trim"
    assert run["backend"] == "tool"
    assert run["command"] == ["tool", "--flag"]
    manifest = json.loads(
        (tmp_path / "microsuite-results.json").read_text(encoding="utf-8")
    )
    assert manifest["schema_version"] == "microsuite-results.v1"
    assert manifest["producer"]["name"] == "microsuite"
    assert manifest["executions"][0]["task"] == "trim"


def test_run_command_writes_results_manifest_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "", ""),
    )

    run_command(
        ["qiime", "diversity-lib", "shannon-entropy"],
        "diversity failed",
        run_dir=tmp_path,
        log=CommandLog(
            task="diversity_calc",
            backend="qiime2",
            inputs={"table": "table.qza"},
            outputs={"vector": "alpha.qza"},
            params={"metric": "shannon"},
        ),
    )

    manifest = json.loads(
        (tmp_path / "microsuite-results.json").read_text(encoding="utf-8")
    )
    assert manifest["run_id"] == tmp_path.name
    assert manifest["executions"][0]["inputs"] == {"table": "table.qza"}
    assert manifest["executions"][0]["params"] == {"metric": "shannon"}
    assert manifest["artifacts"] == [
        {
            "backend": "qiime2",
            "format": "qza",
            "id": "diversity_calc|qiime2|vector|alpha.qza",
            "kind": "alpha_diversity",
            "label": "vector",
            "path": "alpha.qza",
            "task": "diversity_calc",
        }
    ]


def test_run_command_appends_results_manifest_without_duplicate_artifacts(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "", ""),
    )

    log = CommandLog(
        task="denoise",
        backend="qiime2-dada2",
        outputs={"table": "table.qza"},
    )
    run_command(["qiime", "dada2"], "failed", run_dir=tmp_path, log=log)
    run_command(["qiime", "dada2"], "failed", run_dir=tmp_path, log=log)

    manifest = json.loads(
        (tmp_path / "microsuite-results.json").read_text(encoding="utf-8")
    )
    assert len(manifest["executions"]) == 2
    assert len(manifest["artifacts"]) == 1
    assert manifest["artifacts"][0]["kind"] == "feature_table"


def test_run_command_failure_logs_and_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 2, "", "bad\n"),
    )

    with pytest.raises(MicrobiomeSuiteError, match="bad"):
        run_command(["tool"], "fallback", run_dir=tmp_path)

    assert not (tmp_path / "microsuite-results.json").exists()
    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert events[-1]["event"] == "command_end"
    assert events[-1]["exit_code"] == 2


def test_run_command_timeout_logs_and_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def timeout_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        raise subprocess.TimeoutExpired(command, timeout=1)

    monkeypatch.setattr("subprocess.run", timeout_run)

    with pytest.raises(MicrobiomeSuiteError, match="timed out"):
        run_command(["tool"], "fallback", run_dir=tmp_path, timeout=1)

    events = [
        json.loads(line)
        for line in (tmp_path / "events.jsonl").read_text(encoding="utf-8").splitlines()
    ]
    assert events[-1]["event"] == "command_timeout"
    assert (tmp_path / "stdout.log").read_text(encoding="utf-8") == ""
    assert (tmp_path / "stderr.log").read_text(encoding="utf-8") == ""
    run = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert run["status"] == "timeout"
    assert run["timeout"] == 1


def test_resolve_threads_auto_and_bounds(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("os.cpu_count", lambda: 8)

    assert resolve_threads("auto") == 7
    assert resolve_threads("4") == 4
    assert resolve_threads(20) == 20

    with pytest.raises(MicrobiomeSuiteError, match="threads"):
        resolve_threads("none")
