from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.methods.table_io import export_table, normalize_table

app = typer.Typer(help="Transform and export feature/profile tables as TSV.", no_args_is_help=True)


@app.command("export")
def export_cmd(
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output matrix TSV.")],
    layer: Annotated[
        str | None, typer.Option("--layer", help="Export this layer instead of X.")
    ] = None,
    metadata: Annotated[
        Path | None, typer.Option("--metadata", help="Also write obs metadata to this TSV.")
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    export_table(table=table, output=output, layer=layer, metadata=metadata, force=force)


@app.command("normalize")
def normalize_cmd(
    method: Annotated[
        str, typer.Option("--method", help="relative, total-sum, clr, or prevalence-filter.")
    ],
    input_path: Annotated[Path, typer.Option("--input", help="Input matrix TSV.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output matrix TSV.")],
    target_sum: Annotated[float, typer.Option("--target-sum")] = 1_000_000.0,
    pseudocount: Annotated[float, typer.Option("--pseudocount")] = 1.0,
    min_prevalence: Annotated[float, typer.Option("--min-prevalence")] = 0.1,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    normalize_table(
        method=method,
        input_path=input_path,
        output=output,
        target_sum=target_sum,
        pseudocount=pseudocount,
        min_prevalence=min_prevalence,
        force=force,
    )
