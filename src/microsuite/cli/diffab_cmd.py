from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.diffab.ancombc import run_ancombc
from microsuite.io.h5ad import read_h5ad

app = typer.Typer(help="Differential abundance commands.", no_args_is_help=True)


@app.command("ancombc")
def ancombc(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    group: Annotated[str, typer.Option("--group", help="obs column defining groups.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    run_ancombc(adata, group=group, output=prepare_output(output, force=force))
