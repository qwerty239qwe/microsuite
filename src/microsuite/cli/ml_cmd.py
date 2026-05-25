from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.methods.ml_longitudinal import longitudinal, ml_classify

app = typer.Typer(help="Machine learning and longitudinal analysis.", no_args_is_help=True)


@app.command("classify")
def classify_cmd(
    backend: Annotated[str, typer.Option("--backend", help="randomforest or xgboost.")],
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad feature table.")],
    target: Annotated[str, typer.Option("--target", help="Sample metadata target column.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Prediction TSV.")],
    importance_output: Annotated[
        Path | None, typer.Option("--importance-output", help="Feature-importance TSV.")
    ] = None,
    test_fraction: Annotated[float, typer.Option("--test-fraction", min=0.0, max=1.0)] = 0.25,
    n_estimators: Annotated[int, typer.Option("--n-estimators", min=1)] = 100,
    seed: Annotated[int, typer.Option("--seed")] = 0,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
) -> None:
    ml_classify(
        backend=backend,
        table=table,
        target=target,
        output=output,
        importance_output=importance_output,
        test_fraction=test_fraction,
        n_estimators=n_estimators,
        seed=seed,
        force=force,
    )


@app.command("longitudinal")
def longitudinal_cmd(
    backend: Annotated[str, typer.Option("--backend", help="native-time-series.")],
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad feature table.")],
    subject: Annotated[str, typer.Option("--subject", help="Subject/sample unit column.")],
    time: Annotated[str, typer.Option("--time", help="Numeric time column.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    group: Annotated[str | None, typer.Option("--group", help="Optional grouping column.")] = None,
    level: Annotated[str | None, typer.Option("--level", help="Optional taxonomy level.")] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    longitudinal(
        backend=backend,
        table=table,
        subject=subject,
        time=time,
        output=output,
        group=group,
        level=level,
        force=force,
    )
