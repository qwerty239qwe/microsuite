from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.methods.network import network

app = typer.Typer(help="Network inference.", no_args_is_help=True)


@app.command("infer")
def infer_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Network backend.")],
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad feature table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output edge-list TSV.")],
    method: Annotated[str, typer.Option("--method", help="pearson or spearman.")] = "pearson",
    transform: Annotated[
        str, typer.Option("--transform", help="counts, relative, or clr.")
    ] = "relative",
    min_abs_weight: Annotated[
        float,
        typer.Option("--min-abs-weight", min=0.0, max=1.0, help="Minimum absolute edge weight."),
    ] = 0.3,
    min_prevalence: Annotated[
        float,
        typer.Option("--min-prevalence", min=0.0, max=1.0, help="Minimum feature prevalence."),
    ] = 0.1,
    top_n: Annotated[int | None, typer.Option("--top-n", min=1)] = None,
    pseudocount: Annotated[float, typer.Option("--pseudocount", min=0.0)] = 1.0,
    spieceasi_method: Annotated[
        str, typer.Option("--spieceasi-method", help="SPIEC-EASI method: mb or glasso.")
    ] = "mb",
    lambda_min_ratio: Annotated[float, typer.Option("--lambda-min-ratio", min=0.0, max=1.0)] = 0.01,
    nlambda: Annotated[int, typer.Option("--nlambda", min=1)] = 20,
    sensitive: Annotated[bool, typer.Option("--sensitive", help="FlashWeave sensitive mode.")] = (
        False
    ),
    heterogeneous: Annotated[
        bool, typer.Option("--heterogeneous", help="FlashWeave heterogeneous mode.")
    ] = False,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
) -> None:
    network(
        backend=backend,
        table=table,
        output=output,
        method=method,
        transform=transform,
        min_abs_weight=min_abs_weight,
        min_prevalence=min_prevalence,
        top_n=top_n,
        pseudocount=pseudocount,
        spieceasi_method=spieceasi_method,
        lambda_min_ratio=lambda_min_ratio,
        nlambda=nlambda,
        sensitive=sensitive,
        heterogeneous=heterogeneous,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )
