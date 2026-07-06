from __future__ import annotations

import typer
from rich.console import Console

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli import (
    data_cmd,
    diffab_cmd,
    diversity_cmd,
    example_cmd,
    import_cmd,
    method_cmd,
    ml_cmd,
    network_cmd,
    ordination_cmd,
    qiime_cmd,
    refdb_cmd,
    viz_cmd,
    workflow_cmd,
)

console = Console(stderr=True)


app = typer.Typer(
    name="microsuite",
    help="Unified microbiome feature-table analysis CLI.",
    no_args_is_help=True,
)


def _install_groups() -> None:
    app.add_typer(import_cmd.app, name="import")
    app.add_typer(diversity_cmd.app, name="diversity")
    app.add_typer(ordination_cmd.app, name="ordination")
    app.add_typer(diffab_cmd.app, name="diffab")
    app.add_typer(viz_cmd.app, name="viz")
    app.add_typer(data_cmd.app, name="data")
    app.add_typer(example_cmd.app, name="example")
    app.add_typer(qiime_cmd.app, name="qiime")
    app.add_typer(refdb_cmd.app, name="refdb")
    app.add_typer(workflow_cmd.app, name="workflow")
    app.add_typer(network_cmd.app, name="network")
    app.add_typer(ml_cmd.app, name="ml")
    app.add_typer(method_cmd.app)


_install_groups()


def main() -> None:
    try:
        app()
    except MicrobiomeSuiteError as exc:
        console.print(f"[red]Error:[/red] {exc}")
        raise typer.Exit(1) from exc


if __name__ == "__main__":
    main()
