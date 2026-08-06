from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

DEFAULT_DADA2_IMAGE = "ghcr.io/qwerty239qwe/microsuite/r-dada2:latest"
_DADA2_IMAGE_ENV = "MICROSUITE_R_DADA2_IMAGE"
DEFAULT_ECOLOGY_IMAGE = "ghcr.io/qwerty239qwe/microsuite/r-ecology:latest"
DEFAULT_VEGAN_IMAGE = DEFAULT_ECOLOGY_IMAGE
_ECOLOGY_IMAGE_ENV = "MICROSUITE_R_ECOLOGY_IMAGE"
_VEGAN_IMAGE_ENV = "MICROSUITE_R_VEGAN_IMAGE"
DEFAULT_DIFFAB_IMAGE_PREFIX = "ghcr.io/qwerty239qwe/microsuite/r-diffab-"
_DIFFAB_IMAGE_ENV_PREFIX = "MICROSUITE_R_DIFFAB_"
DEFAULT_BATCH_IMAGE_PREFIX = "ghcr.io/qwerty239qwe/microsuite/"
_BATCH_IMAGE_ENV_PREFIX = "MICROSUITE_R_BATCH_"


@dataclass(frozen=True)
class Mount:
    host: Path
    container: str
    mode: str = "rw"


@dataclass(frozen=True)
class EngineProbe:
    available: bool
    responsive: bool
    executable: str | None = None
    version: str | None = None
    error: str | None = None


@dataclass(frozen=True)
class ImageProbe:
    available: bool
    image: str
    digest: str | None = None
    error: str | None = None


@dataclass(frozen=True)
class BindMountProbe:
    writable: bool
    host_path: Path
    image: str
    error: str | None = None


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


def probe_engine(engine: str = "docker", *, timeout: int = 10) -> EngineProbe:
    """Probe an engine executable and daemon without raising user-facing errors."""
    executable = shutil.which(engine)
    if executable is None:
        return EngineProbe(
            available=False,
            responsive=False,
            error=f"'{engine}' was not found on PATH.",
        )
    try:
        result = subprocess.run(
            [executable, "info", "--format", "{{.ServerVersion}}"],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return EngineProbe(True, False, executable=executable, error=str(exc))
    if result.returncode != 0:
        error = result.stderr.strip() or result.stdout.strip() or "engine daemon did not respond"
        return EngineProbe(True, False, executable=executable, error=error)
    return EngineProbe(
        True,
        True,
        executable=executable,
        version=result.stdout.strip() or None,
    )


def probe_image(engine: str, image: str, *, timeout: int = 10) -> ImageProbe:
    """Check whether an image is locally inspectable without pulling it."""
    executable = shutil.which(engine)
    if executable is None:
        return ImageProbe(False, image, error=f"'{engine}' was not found on PATH.")
    try:
        result = subprocess.run(
            [executable, "image", "inspect", "--format", "{{.Id}}", image],
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return ImageProbe(False, image, error=str(exc))
    if result.returncode != 0:
        error = result.stderr.strip() or result.stdout.strip() or "image is unavailable"
        return ImageProbe(False, image, error=error)
    digest = result.stdout.strip() or None
    return ImageProbe(True, image, digest=digest)


def probe_bind_mount(
    engine: str,
    image: str,
    host_dir: Path,
    *,
    timeout: int = 20,
) -> BindMountProbe:
    """Verify that a container can write through a temporary bind mount."""
    executable = shutil.which(engine)
    host_dir = Path(host_dir).resolve()
    if executable is None:
        return BindMountProbe(False, host_dir, image, f"'{engine}' was not found on PATH.")
    try:
        with tempfile.TemporaryDirectory(prefix=".microsuite-bind-", dir=host_dir) as probe_dir:
            command = [executable, "run", "--rm"]
            user = host_user_spec()
            if user is not None:
                command.extend(["--user", user])
            command.extend(
                [
                    "-v",
                    f"{probe_dir}:/microsuite-probe",
                    image,
                    "sh",
                    "-c",
                    "printf ok > /microsuite-probe/write-test",
                ]
            )
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout,
                check=False,
            )
            marker = Path(probe_dir) / "write-test"
            if result.returncode == 0 and marker.read_text(encoding="utf-8") == "ok":
                return BindMountProbe(True, host_dir, image)
            error = result.stderr.strip() or result.stdout.strip() or "bind-mount write failed"
            return BindMountProbe(False, host_dir, image, error)
    except (OSError, subprocess.SubprocessError) as exc:
        return BindMountProbe(False, host_dir, image, str(exc))


def resolve_dada2_image(override: str | None) -> str:
    if override:
        return override
    env = os.environ.get(_DADA2_IMAGE_ENV)
    if env:
        return env
    return DEFAULT_DADA2_IMAGE


def resolve_ecology_image(override: str | None) -> str:
    """Resolve the vegan ecology image: explicit override, env, then GHCR default."""
    if override:
        return override
    env = os.environ.get(_ECOLOGY_IMAGE_ENV) or os.environ.get(_VEGAN_IMAGE_ENV)
    if env:
        return env
    return DEFAULT_ECOLOGY_IMAGE


def resolve_vegan_image(override: str | None) -> str:
    """Compatibility alias for callers that name the backend rather than image."""
    return resolve_ecology_image(override)


def resolve_diffab_image(backend: str, override: str | None) -> str:
    """Resolve the per-backend r-diffab image: override, then env, then default."""
    if override:
        return override
    env = os.environ.get(f"{_DIFFAB_IMAGE_ENV_PREFIX}{backend.upper()}_IMAGE")
    if env:
        return env
    return f"{DEFAULT_DIFFAB_IMAGE_PREFIX}{backend}:latest"


def resolve_batch_image(backend: str, override: str | None, *, image: str | None = None) -> str:
    """Resolve the per-backend r-batch image: override, then env, then default.

    ``backend`` names the env-var lookup (e.g. ``combat-seq`` ->
    ``MICROSUITE_R_BATCH_COMBAT_SEQ_IMAGE``); ``image`` is the actual image
    basename to interpolate into the default GHCR path, since the two are
    not always the same string (``combat-seq`` builds as
    ``r-batch-combatseq``). Defaults to ``backend`` when omitted, for
    callers that don't distinguish the two.
    """
    if override:
        return override
    env_name = f"{_BATCH_IMAGE_ENV_PREFIX}{backend.upper().replace('-', '_')}_IMAGE"
    env = os.environ.get(env_name)
    if env:
        return env
    basename = image if image is not None else backend
    return f"{DEFAULT_BATCH_IMAGE_PREFIX}{basename}:latest"


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
