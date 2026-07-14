from __future__ import annotations

import json
import os
import shlex
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata.stage import active_stage
from microsuite.runtime.results import write_results_manifest


def _note_active_subprocess(
    command: list[str], *, exit_code: int | None, duration_sec: float, status: str
) -> None:
    """Contribute this subprocess to the active stage, if any (independent of run_dir)."""
    record = active_stage()
    if record is not None:
        record.note_subprocess(
            command, exit_code=exit_code, duration_sec=duration_sec, status=status
        )


@dataclass(frozen=True)
class CommandLog:
    task: str | None = None
    backend: str | None = None
    inputs: Any | None = None
    outputs: Any | None = None
    params: dict[str, Any] | None = None


def resolve_threads(threads: int | str, *, reserved: int = 1) -> int:
    cpus = os.cpu_count() or 1
    if isinstance(threads, str):
        if threads == "auto":
            return max(1, cpus - reserved)
        try:
            parsed = int(threads)
        except ValueError as exc:
            raise MicrobiomeSuiteError("--threads must be a positive integer or 'auto'.") from exc
    else:
        parsed = threads
    if parsed < 1:
        raise MicrobiomeSuiteError("--threads must be at least 1.")
    return parsed


def run_command(
    command: list[str],
    failure_message: str,
    *,
    run_dir: Path | None = None,
    log: CommandLog | None = None,
    timeout: float | None = None,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
) -> subprocess.CompletedProcess[str]:
    started = time.time()
    log = log or CommandLog()
    if run_dir is not None:
        run_dir.mkdir(parents=True, exist_ok=True)
        _write_command(run_dir, command)
        _append_event(run_dir, "command_start", command=command, log=log, started_at=started)
        _write_run(run_dir, command=command, log=log, started_at=started)
    try:
        run_kwargs: dict[str, Any] = {"check": False, "text": True, "capture_output": True}
        if timeout is not None:
            run_kwargs["timeout"] = timeout
        if cwd is not None:
            run_kwargs["cwd"] = cwd
        if env is not None:
            run_kwargs["env"] = env
        result = subprocess.run(command, **run_kwargs)
    except subprocess.TimeoutExpired as exc:
        duration_sec = time.time() - started
        _note_active_subprocess(
            command, exit_code=None, duration_sec=duration_sec, status="timed_out"
        )
        if run_dir is not None:
            (run_dir / "stdout.log").write_text(_timeout_stream(exc.stdout), encoding="utf-8")
            (run_dir / "stderr.log").write_text(_timeout_stream(exc.stderr), encoding="utf-8")
            _append_event(
                run_dir,
                "command_timeout",
                command=command,
                log=log,
                started_at=started,
                duration_sec=duration_sec,
                timeout=timeout,
            )
            _write_run(
                run_dir,
                command=command,
                log=log,
                started_at=started,
                duration_sec=duration_sec,
                status="timeout",
                timeout=timeout,
            )
        raise MicrobiomeSuiteError(
            f"Command timed out after {timeout} seconds: {command[0]}"
        ) from exc
    except OSError as exc:
        duration_sec = time.time() - started
        _note_active_subprocess(
            command, exit_code=None, duration_sec=duration_sec, status="launch_failed"
        )
        raise MicrobiomeSuiteError(f"Failed to launch command: {command[0]}: {exc}") from exc

    duration_sec = time.time() - started
    _note_active_subprocess(
        command,
        exit_code=result.returncode,
        duration_sec=duration_sec,
        status="completed" if result.returncode == 0 else "failed",
    )
    if run_dir is not None:
        (run_dir / "stdout.log").write_text(result.stdout or "", encoding="utf-8")
        (run_dir / "stderr.log").write_text(result.stderr or "", encoding="utf-8")
        _append_event(
            run_dir,
            "command_end",
            command=command,
            log=log,
            started_at=started,
            duration_sec=duration_sec,
            exit_code=result.returncode,
        )
        _write_run(
            run_dir,
            command=command,
            log=log,
            started_at=started,
            duration_sec=duration_sec,
            exit_code=result.returncode,
        )
        if result.returncode == 0:
            write_results_manifest(
                run_dir,
                command=command,
                log=_log_payload(log),
                started_at=started,
                duration_sec=duration_sec,
                exit_code=result.returncode,
            )
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or failure_message
        raise MicrobiomeSuiteError(message)
    return result


def _write_command(run_dir: Path, command: list[str]) -> None:
    command_text = " ".join(shlex.quote(part) for part in command)
    (run_dir / "command.txt").write_text(f"{command_text}\n", encoding="utf-8")


def _append_event(
    run_dir: Path,
    event: str,
    *,
    command: list[str],
    log: CommandLog,
    started_at: float,
    **extra: Any,
) -> None:
    payload = {
        "event": event,
        "time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "started_at": started_at,
        "command": command,
        **_log_payload(log),
        **extra,
    }
    with (run_dir / "events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(payload, sort_keys=True) + "\n")


def _write_run(
    run_dir: Path,
    *,
    command: list[str],
    log: CommandLog,
    started_at: float,
    duration_sec: float | None = None,
    exit_code: int | None = None,
    status: str | None = None,
    timeout: float | None = None,
) -> None:
    payload = {
        "started_at": started_at,
        "command": command,
        **_log_payload(log),
    }
    if duration_sec is not None:
        payload["duration_sec"] = duration_sec
    if exit_code is not None:
        payload["exit_code"] = exit_code
    if status is not None:
        payload["status"] = status
    if timeout is not None:
        payload["timeout"] = timeout
    (run_dir / "run.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8"
    )


def _log_payload(log: CommandLog) -> dict[str, Any]:
    return {key: value for key, value in asdict(log).items() if value is not None}


def _timeout_stream(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode(errors="replace")
    return value
