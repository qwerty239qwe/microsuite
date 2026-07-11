from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.io.h5ad import read_h5ad
from microsuite.viz.barplot import taxonomy_barplot
from microsuite.viz.metadata import plot_taxa_by_group

app = typer.Typer(help="Visualization commands.", no_args_is_help=True)


@app.command("barplot")
def barplot(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    level: Annotated[str, typer.Option("--level", help="Taxonomy level.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output PNG.")],
    top_n: Annotated[int, typer.Option("--top-n", min=1)] = 20,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    taxonomy_barplot(adata, level=level, output=prepare_output(output, force=force), top_n=top_n)


@app.command("taxa-by-group")
def taxa_by_group(
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
    level: Annotated[str, typer.Option("--level", help="Taxonomy level.")],
    group_by: Annotated[str, typer.Option("--group-by", help="Metadata (obs) column to group by.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output PNG.")],
    top_n: Annotated[int, typer.Option("--top-n", min=1)] = 15,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    plot_taxa_by_group(
        adata,
        level=level,
        group_by=group_by,
        output=prepare_output(output, force=force),
        top_n=top_n,
    )
