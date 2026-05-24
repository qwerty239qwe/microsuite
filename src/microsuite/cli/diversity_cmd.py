from __future__ import annotations

from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity.alpha import alpha_diversity
from microsuite.diversity.beta import beta_diversity
from microsuite.diversity.ecology import (
    beta_significance,
    beta_turnover,
    constrained_ordination,
    gamma_diversity,
    mantel_test,
    taxa_turnover,
)
from microsuite.io.h5ad import read_h5ad

app = typer.Typer(help="Diversity and ecological statistics.", no_args_is_help=True)


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


@app.command("beta-significance")
def beta_significance_cmd(
    distance_matrix: Annotated[Path, typer.Argument(help="Input square distance matrix TSV.")],
    metadata: Annotated[Path, typer.Option("--metadata", "-m", help="Sample metadata TSV.")],
    column: Annotated[str, typer.Option("--column", help="Metadata grouping column.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    method: Annotated[str, typer.Option("--method", help="permanova or anosim.")] = "permanova",
    permutations: Annotated[int, typer.Option("--permutations", min=0)] = 999,
    seed: Annotated[int, typer.Option("--seed")] = 0,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    result = beta_significance(
        _read_distance_matrix(distance_matrix),
        _read_metadata(metadata),
        column=column,
        method=method,
        permutations=permutations,
        seed=seed,
    )
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


@app.command("mantel")
def mantel_cmd(
    matrix_a: Annotated[Path, typer.Argument(help="First square distance matrix TSV.")],
    matrix_b: Annotated[Path, typer.Argument(help="Second square distance matrix TSV.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    method: Annotated[str, typer.Option("--method", help="pearson or spearman.")] = "pearson",
    permutations: Annotated[int, typer.Option("--permutations", min=0)] = 999,
    seed: Annotated[int, typer.Option("--seed")] = 0,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    result = mantel_test(
        _read_distance_matrix(matrix_a),
        _read_distance_matrix(matrix_b),
        method=method,
        permutations=permutations,
        seed=seed,
    )
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


@app.command("gamma")
def gamma_cmd(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    group: Annotated[str, typer.Option("--group", help="Sample metadata group column.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    metric: Annotated[str, typer.Option("--metric", help="Alpha metric for pooled groups.")] = (
        "observed_features"
    ),
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    result = gamma_diversity(read_h5ad(ensure_input(table)), group=group, metric=metric)
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


@app.command("beta-turnover")
def beta_turnover_cmd(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    level: Annotated[str | None, typer.Option("--level", help="Optional taxonomy level.")] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    result = beta_turnover(read_h5ad(ensure_input(table)), level=level)
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


@app.command("taxa-turnover")
def taxa_turnover_cmd(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    group: Annotated[str, typer.Option("--group", help="Sample metadata group column.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    level: Annotated[str | None, typer.Option("--level", help="Optional taxonomy level.")] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    result = taxa_turnover(read_h5ad(ensure_input(table)), group=group, level=level)
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


@app.command("constrained-ordination")
def constrained_ordination_cmd(
    table: Annotated[Path, typer.Argument(help="Input .h5ad table.")],
    constraint: Annotated[
        list[str],
        typer.Option("--constraint", help="Sample metadata constraint. Repeat as needed."),
    ],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    method: Annotated[str, typer.Option("--method", help="rda, db-rda, or cca.")] = "rda",
    dimensions: Annotated[int, typer.Option("--dimensions", min=1)] = 2,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    result = constrained_ordination(
        read_h5ad(ensure_input(table)),
        constraints=constraint,
        method=method,
        dimensions=dimensions,
    )
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


def _read_distance_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(ensure_input(path), sep="\t", index_col=0)


def _read_metadata(path: Path) -> pd.DataFrame:
    return pd.read_csv(ensure_input(path), sep="\t", index_col=0)
