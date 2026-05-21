from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output


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


def run_qiime(command: list[str], failure_message: str) -> None:
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or failure_message
        raise MicrobiomeSuiteError(message)
