from __future__ import annotations

import os
import shutil
import sys
import tempfile
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from platformdirs import user_cache_path

from microsuite.runtime import container
from microsuite.system.capabilities import require_capabilities
from microsuite.system.version import package_version

CheckStatus = Literal["pass", "warn", "fail", "skip"]


@dataclass(frozen=True)
class DoctorCheck:
    id: str
    status: CheckStatus
    message: str
    details: dict[str, Any] | None = None
    remediation: str | None = None


@dataclass(frozen=True)
class DoctorReport:
    checks: tuple[DoctorCheck, ...]

    @property
    def exit_code(self) -> int:
        return 1 if any(check.status == "fail" for check in self.checks) else 0

    def to_dict(self) -> dict[str, Any]:
        status = "fail" if self.exit_code else "pass"
        if status == "pass" and any(check.status == "warn" for check in self.checks):
            status = "warn"
        return {
            "schema_version": "microsuite-doctor.v1",
            "producer": {"name": "microsuite", "version": package_version()},
            "status": status,
            "checks": [asdict(check) for check in self.checks],
        }


def available_memory_bytes() -> int | None:
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def _writable_path_check(identifier: str, path: Path) -> DoctorCheck:
    try:
        path.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(prefix=".microsuite-write-", dir=path):
            pass
    except OSError as exc:
        return DoctorCheck(
            identifier,
            "fail",
            f"Cannot write to {path}",
            details={"path": str(path), "error": str(exc)},
            remediation="Choose a writable path or correct its ownership and permissions.",
        )
    return DoctorCheck(identifier, "pass", f"Writable: {path}", details={"path": str(path)})


def _memory_check(minimum_gb: float) -> DoctorCheck:
    available = available_memory_bytes()
    if available is None:
        return DoctorCheck("host.memory", "warn", "Available memory could not be determined.")
    available_gb = available / 1024**3
    status: CheckStatus = "pass" if available_gb >= minimum_gb else "warn"
    return DoctorCheck(
        "host.memory",
        status,
        f"{available_gb:.1f} GiB memory available.",
        details={"available_bytes": available, "minimum_gb": minimum_gb},
        remediation=(
            None if status == "pass" else "Increase available memory for memory-intensive stages."
        ),
    )


def _disk_check(path: Path, minimum_gb: float) -> DoctorCheck:
    try:
        free = shutil.disk_usage(path).free
    except OSError as exc:
        return DoctorCheck(
            "host.disk",
            "warn",
            f"Free disk space could not be determined for {path}.",
            details={"error": str(exc)},
        )
    free_gb = free / 1024**3
    status: CheckStatus = "pass" if free_gb >= minimum_gb else "warn"
    return DoctorCheck(
        "host.disk",
        status,
        f"{free_gb:.1f} GiB free at {path}.",
        details={"free_bytes": free, "minimum_gb": minimum_gb},
    )


def run_doctor(
    *,
    output_dir: Path | None = None,
    cache_dir: Path | None = None,
    engine: str = "docker",
    required_capabilities: Iterable[str] = (),
    executables: Iterable[str] = (),
    images: Iterable[str] = (),
    require_container: bool = False,
    minimum_memory_gb: float = 4.0,
    minimum_disk_gb: float = 5.0,
) -> DoctorReport:
    output = Path.cwd() if output_dir is None else Path(output_dir)
    cache = Path(user_cache_path("microsuite")) if cache_dir is None else Path(cache_dir)
    required_capabilities = tuple(required_capabilities)
    executables = tuple(executables)
    images = tuple(images)
    checks: list[DoctorCheck] = [
        DoctorCheck(
            "host.python",
            "pass",
            f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            details={"executable": sys.executable},
        ),
        DoctorCheck(
            "host.user",
            "pass",
            "Host user identity detected.",
            details={
                "uid": getattr(os, "getuid", lambda: None)(),
                "gid": getattr(os, "getgid", lambda: None)(),
            },
        ),
        _writable_path_check("path.output", output),
        _writable_path_check("path.cache", cache),
        _memory_check(minimum_memory_gb),
        _disk_check(output, minimum_disk_gb),
    ]

    engine_probe = container.probe_engine(engine)
    if engine_probe.responsive:
        checks.append(
            DoctorCheck(
                "container.engine",
                "pass",
                f"{engine} daemon is available.",
                details={"executable": engine_probe.executable, "version": engine_probe.version},
            )
        )
    else:
        status: CheckStatus = "fail" if require_container or images else "warn"
        checks.append(
            DoctorCheck(
                "container.engine",
                status,
                engine_probe.error or f"{engine} is unavailable.",
                remediation=(
                    f"Install {engine} and ensure the current user can access its daemon."
                ),
            )
        )

    for executable in executables:
        resolved = shutil.which(executable)
        message = (
            f"Executable found: {resolved}" if resolved else f"Executable not found: {executable}"
        )
        checks.append(
            DoctorCheck(
                f"executable.{executable}",
                "pass" if resolved else "fail",
                message,
                details={"path": resolved},
                remediation=None if resolved else f"Install {executable} and add it to PATH.",
            )
        )

    available, missing = require_capabilities(required_capabilities)
    for capability in available:
        checks.append(DoctorCheck(f"capability.{capability}", "pass", "Capability is available."))
    for capability in missing:
        checks.append(
            DoctorCheck(
                f"capability.{capability}",
                "fail",
                "Required capability is unavailable.",
                remediation="Update MicroSuite or select a supported workflow implementation.",
            )
        )

    if engine_probe.responsive:
        for image in images:
            probe = container.probe_image(engine, image)
            message = (
                f"Image available: {image}" if probe.available else f"Image unavailable: {image}"
            )
            remediation = None if probe.available else f"Pull the image with: {engine} pull {image}"
            checks.append(
                DoctorCheck(
                    f"container.image.{image}",
                    "pass" if probe.available else "fail",
                    message,
                    details={"digest": probe.digest, "error": probe.error},
                    remediation=remediation,
                )
            )
            if probe.available:
                bind_probe = container.probe_bind_mount(engine, image, output)
                checks.append(
                    DoctorCheck(
                        f"container.bind_mount.{image}",
                        "pass" if bind_probe.writable else "fail",
                        (
                            f"Container can write to {output}."
                            if bind_probe.writable
                            else f"Container cannot write to {output}."
                        ),
                        details={"path": str(output), "error": bind_probe.error},
                        remediation=(
                            None
                            if bind_probe.writable
                            else "Check bind-mount sharing and host directory permissions."
                        ),
                    )
                )

    return DoctorReport(tuple(checks))
