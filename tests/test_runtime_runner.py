from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command


def test_run_command_writes_structured_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls: list[list[str]] = []

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
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
    manifest = json.loads((tmp_path / "microsuite-results.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "microsuite-results.v1"
    assert manifest["producer"]["name"] == "microsuite"
    assert manifest["executions"][0]["task"] == "trim"


def test_run_command_duration_uses_monotonic_clock_when_wall_clock_moves_backward(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    wall_times = iter([1000.0])
    monotonic_times = iter([50.0, 55.0])
    monkeypatch.setattr("microsuite.runtime.runner.time.time", lambda: next(wall_times))
    monkeypatch.setattr("microsuite.runtime.runner.time.monotonic", lambda: next(monotonic_times))
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "", ""),
    )

    run_command(["tool"], "failed", run_dir=tmp_path)

    run = json.loads((tmp_path / "run.json").read_text(encoding="utf-8"))
    assert run["duration_sec"] == 5.0


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

    manifest = json.loads((tmp_path / "microsuite-results.json").read_text(encoding="utf-8"))
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

    manifest = json.loads((tmp_path / "microsuite-results.json").read_text(encoding="utf-8"))
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


def test_run_command_survives_output_that_is_not_valid_utf8(tmp_path: Path) -> None:
    # External tools are not obliged to emit valid UTF-8. mothur's align.seqs
    # draws a progress bar containing raw bytes, and bare text=True decoded
    # strictly, so the whole run died with UnicodeDecodeError instead of
    # executing the tool. Undecodable bytes must degrade, not abort.
    result = run_command(
        [
            sys.executable,
            "-c",
            r"import sys; sys.stdout.buffer.write(b'start\xce\xff\xfeend')",
        ],
        "should not fail",
        run_dir=tmp_path / "run",
    )

    assert "start" in result.stdout
    assert "end" in result.stdout
    # The captured log must also be writable, which it would not be if the
    # replacement characters were absent and raw bytes had survived.
    assert "start" in (tmp_path / "run" / "stdout.log").read_text(encoding="utf-8")


def test_run_command_decodes_output_as_utf8_regardless_of_platform_locale(
    tmp_path: Path,
) -> None:
    # text=True without an explicit encoding uses the platform locale, so the
    # same tool bytes decoded as cp1252 on Windows and UTF-8 on Linux.
    result = run_command(
        [sys.executable, "-c", r"import sys; sys.stdout.buffer.write('café·µ'.encode())"],
        "should not fail",
    )

    assert "café·µ" in result.stdout
