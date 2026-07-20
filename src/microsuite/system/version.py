from __future__ import annotations

import platform
import subprocess
from importlib import metadata
from pathlib import Path
from typing import Any


def package_version() -> str:
    """Return the installed distribution version without importing microsuite."""
    try:
        return metadata.version("microsuite")
    except metadata.PackageNotFoundError:
        return "0+unknown"


def source_commit() -> tuple[str, str | None]:
    """Report whether this is a source checkout and its best-effort commit ID."""
    repository = Path(__file__).resolve().parents[3]
    if not (repository / ".git").exists():
        return "installed", None
    try:
        result = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=repository,
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return "editable", None
    commit = result.stdout.strip() if result.returncode == 0 else None
    return "editable", commit or None


def version_info() -> dict[str, Any]:
    source, commit = source_commit()
    return {
        "name": "microsuite",
        "version": package_version(),
        "source": source,
        "commit": commit,
        "python": platform.python_version(),
    }
