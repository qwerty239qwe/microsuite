"""Shared local/docker runner for the R differential-abundance backends.

Both wrappers (``ancombc.py`` and ``r_backends.py``) build an ordered positional
argument list for their ``.R`` script and hand it here. Convention: the last
``Path`` is the output (its parent is bind-mounted read-write); every earlier
``Path`` is an input (its parent is bind-mounted read-only); ``str`` items pass
through verbatim. Docker runs execute the per-backend image as the caller's
UID/GID (writable outputs) and record the resolved image + digest in
``<output-dir>/<backend>_container.json``.
"""

from __future__ import annotations

import json
import shutil
from importlib.resources import files
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.container import (
    PathMapper,
    build_container_command,
    host_user_spec,
    require_engine,
    resolve_diffab_image,
    resolve_image_digest,
)
from microsuite.runtime.runner import CommandLog, run_command


def invoke_r_backend(
    *,
    backend: str,
    positional: list[str | Path],
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
    run_dir: Path | None = None,
    timeout: float | None = None,
    log: CommandLog,
    local_missing_message: str,
) -> None:
    if runtime not in ("local", "docker"):
        raise MicrobiomeSuiteError(
            f"Unsupported --runtime '{runtime}' for {backend}; choose 'local' or 'docker'."
        )

    if runtime == "local":
        rscript = shutil.which("Rscript")
        if rscript is None:
            raise MicrobiomeSuiteError(local_missing_message)
        script = files("microsuite.diffab.r").joinpath(f"{backend}.R")
        command = [rscript, str(script), *[str(arg) for arg in positional]]
        run_command(
            command,
            failure_message=f"{backend} failed.",
            run_dir=run_dir,
            log=log,
            timeout=timeout,
        )
        return

    resolved_image = resolve_diffab_image(backend, image)
    require_engine(engine)
    paths = [arg for arg in positional if isinstance(arg, Path)]
    if not paths:
        raise MicrobiomeSuiteError(f"{backend} requires at least one file argument.")
    output_path = paths[-1]

    mapper = PathMapper()
    mountpoints: dict[Path, str] = {}

    def _mount(host_dir: Path, mode: str) -> None:
        resolved = host_dir.resolve()
        if resolved not in mountpoints:
            mountpoints[resolved] = f"/mnt/d{len(mountpoints)}"
        mapper.add_dir(host_dir, mode, mountpoints[resolved])

    for path in paths[:-1]:
        _mount(path.parent, "ro")
    _mount(output_path.parent, "rw")

    inner = [f"/opt/microsuite/{backend}.R"]
    for arg in positional:
        inner.append(mapper.to_container(arg) if isinstance(arg, Path) else arg)
    command = build_container_command(
        inner, resolved_image, mapper.mounts(), engine=engine, user=host_user_spec()
    )
    run_command(
        command,
        failure_message=f"{backend} failed.",
        run_dir=run_dir,
        log=log,
        timeout=timeout,
    )

    sidecar = output_path.parent / f"{backend}_container.json"
    sidecar.write_text(
        json.dumps(
            {
                "runtime": "docker",
                "engine": engine,
                "image": resolved_image,
                "digest": resolve_image_digest(engine, resolved_image),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
