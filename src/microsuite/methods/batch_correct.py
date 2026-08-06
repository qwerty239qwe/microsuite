from __future__ import annotations

from pathlib import Path

from microsuite._paths import ensure_input, prepare_output
from microsuite.batch.backends import SUPPORTED_BACKENDS
from microsuite.batch.correct import run_batch_correction
from microsuite.io.h5ad import read_h5ad, write_h5ad

__all__ = ["SUPPORTED_BACKENDS", "batch_correct"]


def batch_correct(
    *,
    backend: str,
    table: Path,
    output: Path,
    batch: str,
    covariates: list[str] | None = None,
    target: str | None = None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
) -> None:
    adata = read_h5ad(ensure_input(table))
    corrected = run_batch_correction(
        adata,
        backend=backend,
        batch=batch,
        covariates=covariates,
        target=target,
        run_dir=run_dir,
        timeout=timeout,
        runtime=runtime,
        image=image,
        engine=engine,
    )
    write_h5ad(corrected, prepare_output(output, force=force))
