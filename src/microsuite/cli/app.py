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
    system_cmd,
    table_cmd,
    viz_cmd,
    workflow_cmd,
)

console = Console(stderr=True)


app = typer.Typer(
    name="microsuite",
    help="Unified microbiome feature-table analysis CLI.",
    no_args_is_help=True,
)

_DEBUG_ENABLED = False


@app.callback()
def root(
    debug: bool = typer.Option(False, "--debug", help="Show tracebacks for known errors."),
) -> None:
    global _DEBUG_ENABLED
    _DEBUG_ENABLED = debug


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
    app.add_typer(table_cmd.app, name="table")
    app.add_typer(method_cmd.app)
    app.command("version")(system_cmd.version_command)
    app.command("capabilities")(system_cmd.capabilities_command)
    app.command("doctor")(system_cmd.doctor_command)


_install_groups()


def main() -> None:
    try:
        app()
    except MicrobiomeSuiteError as exc:
        if _DEBUG_ENABLED:
            raise
        console.print(f"[red]Error:[/red] {exc}")
        raise SystemExit(1) from None


if __name__ == "__main__":
    main()
