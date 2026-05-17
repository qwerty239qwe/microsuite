from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity.alpha import alpha_diversity
from microsuite.diversity.beta import beta_diversity
from microsuite.io.h5ad import read_h5ad

app = typer.Typer(help="Diversity metrics.", no_args_is_help=True)


@app.command("alpha")
def alpha(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    metric: Annotated[str, typer.Option("--metric", help="Alpha metric.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    result = alpha_diversity(adata, metric)
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


@app.command("beta")
def beta(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    metric: Annotated[str, typer.Option("--metric", help="Beta metric.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output square TSV matrix.")],
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    result = beta_diversity(adata, metric)
    result.to_csv(prepare_output(output, force=force), sep="\t")
