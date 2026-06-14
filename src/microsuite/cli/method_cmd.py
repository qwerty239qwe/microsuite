from __future__ import annotations

import typer

from microsuite.cli import (
    method_diversity_cmd,
    method_features_cmd,
    method_preprocess_cmd,
    method_qiime_io_cmd,
    method_stats_cmd,
    method_tables_cmd,
    method_taxonomy_cmd,
)
from microsuite.cli._method_api import METHOD_BACKENDS, QIIME2_WRAPPER_METHODS

app = typer.Typer(help="Method-oriented microbiome operations.", no_args_is_help=True)


@app.command("methods")
def methods() -> None:
    for command, backends in METHOD_BACKENDS.items():
        typer.echo(command)
        for backend in backends:
            typer.echo(f"  - {backend}")
    for method, backends in QIIME2_WRAPPER_METHODS.items():
        typer.echo(method)
        for backend in backends:
            typer.echo(f"  - {backend}")


for _module in (
    method_preprocess_cmd,
    method_features_cmd,
    method_taxonomy_cmd,
    method_tables_cmd,
    method_diversity_cmd,
    method_qiime_io_cmd,
    method_stats_cmd,
):
    _module.register(app)
