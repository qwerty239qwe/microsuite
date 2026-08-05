from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.batch.backends import SUPPORTED_BACKENDS
from microsuite.methods.batch_correct import batch_correct

app = typer.Typer(help="Batch effect correction commands.", no_args_is_help=True)


@app.command("correct")
def correct(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output .h5ad table.")],
    batch: Annotated[str, typer.Option("--batch-col", help="obs column holding the batch label.")],
    backend: Annotated[
        str,
        typer.Option("--backend", help=f"One of: {', '.join(SUPPORTED_BACKENDS)}."),
    ] = "mmuphin",
    covariates: Annotated[
        list[str] | None,
        typer.Option("--covariates", help="obs column to preserve (repeatable)."),
    ] = None,
    target: Annotated[
        str | None,
        typer.Option(
            "--target-col",
            help="Outcome column. Required by supervised backends; see "
            "docs/batch_correction.md for the leakage hazard.",
        ),
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
    runtime: Annotated[
        str, typer.Option("--runtime", help="R backend runtime: 'local' Rscript or 'docker'.")
    ] = "local",
    image: Annotated[
        str | None, typer.Option("--image", help="Override the r-batch-<backend> image.")
    ] = None,
) -> None:
    batch_correct(
        backend=backend,
        table=table,
        output=output,
        batch=batch,
        covariates=list(covariates or []),
        target=target,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
        runtime=runtime,
        image=image,
    )
