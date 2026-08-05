"""Compatibility shim: the diffab R runner now lives in ``runtime.r_backend``.

Kept so every existing diffab caller and its tests are unaffected by the move.
New backends should call ``invoke_r_script`` directly.
"""

from __future__ import annotations

from pathlib import Path

from microsuite.runtime.container import resolve_diffab_image
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog


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
    invoke_r_script(
        backend=backend,
        script_package="microsuite.diffab.r",
        script_name=backend,
        resolve_image=resolve_diffab_image,
        positional=positional,
        runtime=runtime,
        image=image,
        engine=engine,
        run_dir=run_dir,
        timeout=timeout,
        log=log,
        local_missing_message=local_missing_message,
    )
