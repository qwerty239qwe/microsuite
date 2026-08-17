from __future__ import annotations

from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity.adonis import adonis2
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
    tree: Annotated[
        Path | None, typer.Option("--tree", help="Newick tree for phylogenetic metrics.")
    ] = None,
    q: Annotated[
        str, typer.Option("--q", help="Comma-separated iNEXT diversity orders.")
    ] = "0,1,2",
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    kwargs = {}
    if metric.lower().replace("-", "_") == "inext":
        kwargs["q"] = _parse_float_csv(q, "--q")
    result = alpha_diversity(
        adata, metric, tree=ensure_input(tree) if tree is not None else None, **kwargs
    )
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
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    column: Annotated[
        str | None,
        typer.Option("--column", help="Metadata grouping column (required by native/anosim2)."),
    ] = None,
    method: Annotated[
        str,
        typer.Option("--method", help="Native: permanova/permdisp/anosim; vegan: adonis2/anosim2."),
    ] = "permanova",
    backend: Annotated[
        str, typer.Option("--backend", help="Significance backend: native or vegan.")
    ] = "native",
    formula: Annotated[
        str | None, typer.Option("--formula", help="R formula right-hand side for vegan.")
    ] = None,
    strata: Annotated[
        str | None,
        typer.Option(
            "--strata",
            help="Vegan permutation block: metadata column or colon-separated interaction.",
        ),
    ] = None,
    permutations: Annotated[int, typer.Option("--permutations", min=0)] = 999,
    seed: Annotated[int, typer.Option("--seed")] = 0,
    runtime: Annotated[
        str, typer.Option("--runtime", help="vegan execution: local Rscript or docker.")
    ] = "local",
    image: Annotated[
        str | None, typer.Option("--image", help="Override the r-ecology container image.")
    ] = None,
    engine: Annotated[
        str, typer.Option("--engine", help="Container engine for --runtime docker.")
    ] = "docker",
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    output_path = prepare_output(output, force=force)
    result = beta_significance(
        _read_distance_matrix(distance_matrix),
        _read_metadata(metadata),
        column=column,
        method=method,
        permutations=permutations,
        seed=seed,
        backend=backend,
        formula=formula,
        strata=strata,
        runtime=runtime,
        image=image,
        engine=engine,
        run_dir=run_dir,
        timeout=timeout,
        sidecar_dir=output_path.parent,
    )
    result.to_csv(output_path, sep="\t", index=False)


@app.command("adonis")
def adonis_cmd(
    distance_matrix: Annotated[Path, typer.Argument(help="Input square distance matrix TSV.")],
    metadata: Annotated[Path, typer.Option("--metadata", "-m", help="Sample metadata TSV.")],
    formula: Annotated[
        str,
        typer.Option(
            "--formula",
            help=(
                "Model formula, e.g. 'dist ~ disease_status + accession"
                " + accession:timepoint + (1 | accession:subject_id)'."
                " A '(1 | group)' term restricts permutation instead of adding coefficients."
            ),
        ),
    ],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    permutations: Annotated[int, typer.Option("--permutations", min=0)] = 999,
    seed: Annotated[int, typer.Option("--seed")] = 0,
    blocks: Annotated[
        str | None,
        typer.Option("--blocks", help="Metadata column samples may never be permuted across."),
    ] = None,
    within: Annotated[
        str,
        typer.Option("--within", help="Shuffle samples inside a group: 'free' or 'none'."),
    ] = "free",
    by: Annotated[
        str,
        typer.Option(
            "--by",
            help=(
                "'terms' for sequential (type I) sums of squares, where order matters; "
                "'margin' to adjust every term for all others so order does not."
            ),
        ),
    ] = "terms",
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    """Multi-term PERMANOVA (vegan adonis2)."""
    result = adonis2(
        _read_distance_matrix(distance_matrix),
        _read_metadata(metadata),
        formula=formula,
        permutations=permutations,
        seed=seed,
        blocks=blocks,
        within=within,
        by=by,
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


def _parse_float_csv(value: str, option_name: str) -> tuple[float, ...]:
    try:
        parsed = tuple(float(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        message = f"{option_name} must be a comma-separated list of numbers."
        raise typer.BadParameter(message) from exc
    if not parsed:
        raise typer.BadParameter(f"{option_name} must contain at least one number.")
    return parsed
