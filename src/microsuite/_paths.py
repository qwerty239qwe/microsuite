from __future__ import annotations

from pathlib import Path

from ._errors import MicrobiomeSuiteError


def ensure_input(path: Path) -> Path:
    if not path.exists():
        raise MicrobiomeSuiteError(f"Input file does not exist: {path}")
    if not path.is_file():
        raise MicrobiomeSuiteError(f"Input path is not a file: {path}")
    return path


def prepare_output(path: Path, *, force: bool) -> Path:
    if path.exists() and not force:
        raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    return path
