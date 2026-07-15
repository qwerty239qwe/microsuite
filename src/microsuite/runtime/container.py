from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

DEFAULT_DADA2_IMAGE = "ghcr.io/qwerty239qwe/microsuite/r-dada2:latest"
_DADA2_IMAGE_ENV = "MICROSUITE_R_DADA2_IMAGE"
DEFAULT_DIFFAB_IMAGE_PREFIX = "ghcr.io/qwerty239qwe/microsuite/r-diffab-"
_DIFFAB_IMAGE_ENV_PREFIX = "MICROSUITE_R_DIFFAB_"


@dataclass(frozen=True)
class Mount:
    host: Path
    container: str
    mode: str = "rw"


def build_container_command(
    inner: list[str],
    image: str,
    mounts: list[Mount],
    *,
    engine: str = "docker",
    user: str | None = None,
) -> list[str]:
    command = [engine, "run", "--rm"]
    if user is not None:
        command.extend(["--user", user])
    for mount in mounts:
        spec = f"{mount.host}:{mount.container}"
        if mount.mode == "ro":
            spec += ":ro"
        command.extend(["-v", spec])
    command.append(image)
    command.extend(inner)
    return command


def host_user_spec() -> str | None:
    """Return the invoking POSIX UID:GID for ownership-safe bind mounts."""
    getuid = getattr(os, "getuid", None)
    getgid = getattr(os, "getgid", None)
    if getuid is None or getgid is None:
        return None
    return f"{getuid()}:{getgid()}"


def require_engine(engine: str = "docker") -> str:
    resolved = shutil.which(engine)
    if resolved is None:
        raise MicrobiomeSuiteError(
            f"The '{engine}' container engine is required for --runtime docker but was "
            f"not found on PATH. Install {engine}, or use --runtime local with the "
            "required backend package installed."
        )
    return resolved


def resolve_dada2_image(override: str | None) -> str:
    if override:
        return override
    env = os.environ.get(_DADA2_IMAGE_ENV)
    if env:
        return env
    return DEFAULT_DADA2_IMAGE


def resolve_diffab_image(backend: str, override: str | None) -> str:
    """Resolve the per-backend r-diffab image: override, then env, then default."""
    if override:
        return override
    env = os.environ.get(f"{_DIFFAB_IMAGE_ENV_PREFIX}{backend.upper()}_IMAGE")
    if env:
        return env
    return f"{DEFAULT_DIFFAB_IMAGE_PREFIX}{backend}:latest"


def resolve_image_digest(engine: str, image: str) -> str | None:
    """Best-effort image digest via `<engine> inspect`; None if unavailable."""
    exe = shutil.which(engine)
    if exe is None:
        return None
    for fmt in ("{{index .RepoDigests 0}}", "{{.Id}}"):
        try:
            result = subprocess.run(
                [exe, "inspect", "--format", fmt, image],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        if result.returncode == 0:
            value = result.stdout.strip()
            if value:
                return value
    return None


class PathMapper:
    """Assign host directories stable container mountpoints and rewrite paths."""

    def __init__(self) -> None:
        self._dirs: dict[Path, Mount] = {}

    def add_dir(self, host_dir: Path, mode: str, container: str) -> None:
        resolved = host_dir.resolve()
        existing = self._dirs.get(resolved)
        if existing is None:
            self._dirs[resolved] = Mount(host=resolved, container=container, mode=mode)
        elif existing.mode == "ro" and mode == "rw":
            self._dirs[resolved] = Mount(host=resolved, container=existing.container, mode="rw")

    def container_dir(self, host_dir: Path) -> str:
        return self._dirs[host_dir.resolve()].container

    def to_container(self, host_path: Path) -> str:
        resolved = host_path.resolve()
        mount = self._dirs.get(resolved.parent)
        if mount is None:
            raise MicrobiomeSuiteError(f"No container mount registered for {resolved.parent}")
        return f"{mount.container}/{resolved.name}"

    def mounts(self) -> list[Mount]:
        return list(self._dirs.values())
