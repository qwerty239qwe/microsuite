from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.runtime.runner import CommandLog, run_command


def require_qiime(task: str) -> str:
    qiime = shutil.which("qiime")
    if qiime is None:
        raise MicrobiomeSuiteError(
            f"{task} requires the external 'qiime' command. "
            "Activate a QIIME 2 environment and rerun this command."
        )
    return qiime


def ensure_inputs(*paths: Path | None) -> None:
    for path in paths:
        if path is not None:
            ensure_input(path)


def prepare_outputs(*paths: Path | None, force: bool) -> None:
    for path in paths:
        if path is not None:
            prepare_output(path, force=force)


def run_qiime(
    command: list[str],
    failure_message: str,
    *,
    run_dir: Path | None = None,
    timeout: float | None = None,
    task: str | None = None,
    backend: str | None = None,
) -> None:
    run_command(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task=task, backend=backend),
    )
