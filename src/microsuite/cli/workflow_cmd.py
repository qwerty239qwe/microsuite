from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.data.moving_pictures import copy_small_fixture
from microsuite.workflows.catalog import WORKFLOWS
from microsuite.workflows.table_summary import run_table_summary

app = typer.Typer(help="Workflow-oriented toolbox commands.", no_args_is_help=True)


@app.command("list")
def list_workflows() -> None:
    for workflow in WORKFLOWS:
        typer.echo(f"{workflow.name}\t{workflow.status}\t{workflow.summary}")


@app.command("show")
def show_workflow(name: Annotated[str, typer.Argument(help="Workflow name.")]) -> None:
    for workflow in WORKFLOWS:
        if workflow.name == name:
            typer.echo(f"name: {workflow.name}")
            typer.echo(f"status: {workflow.status}")
            typer.echo(f"summary: {workflow.summary}")
            typer.echo(f"inputs: {workflow.inputs}")
            typer.echo(f"outputs: {workflow.outputs}")
            return
    raise typer.BadParameter(f"Unknown workflow: {name}")


@app.command("table-summary")
def table_summary(
    output: Annotated[Path, typer.Option("--out", "-o", help="Output run directory.")],
    table: Annotated[Path, typer.Option("--table", help="Input feature table.")],
    metadata: Annotated[Path, typer.Option("--metadata", "-m", help="Sample metadata TSV.")],
    taxonomy: Annotated[Path | None, typer.Option("--taxonomy", "-t")] = None,
    taxonomy_artifact: Annotated[
        Path | None, typer.Option("--taxonomy-artifact", help="QIIME 2 taxonomy artifact.")
    ] = None,
    input_format: Annotated[
        str, typer.Option("--format", help="Input format: tsv or qza.")
    ] = "tsv",
    alpha_metric: Annotated[str, typer.Option("--alpha-metric")] = "shannon",
    beta_metric: Annotated[str, typer.Option("--beta-metric")] = "bray-curtis",
    barplot_level: Annotated[str, typer.Option("--barplot-level")] = "genus",
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
) -> None:
    run_table_summary(
        output=output,
        table=table,
        metadata=metadata,
        taxonomy=taxonomy,
        taxonomy_artifact=taxonomy_artifact,
        input_format=input_format,
        alpha_metric=alpha_metric,
        beta_metric=beta_metric,
        barplot_level=barplot_level,
        force=force,
    )


@app.command("moving-pictures")
def moving_pictures(
    output: Annotated[Path, typer.Option("--out", "-o", help="Output run directory.")],
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
) -> None:
    data_dir = output / "data"
    copy_small_fixture(data_dir, force=force)
    run_table_summary(
        output=output,
        table=data_dir / "table.tsv",
        metadata=data_dir / "metadata.tsv",
        taxonomy=data_dir / "taxonomy.tsv",
        input_format="tsv",
        force=force,
    )
