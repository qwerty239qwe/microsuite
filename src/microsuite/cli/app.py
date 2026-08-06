from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli import (
    batch_cmd,
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
from microsuite.primer import check_fastq_primers, primer_check_fails

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
    app.add_typer(batch_cmd.app, name="batch")
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


@app.command("primer-check")
def primer_check_command(
    inputs: Annotated[
        list[Path],
        typer.Option("--input", help="FASTQ file to inspect. Repeat for multiple files."),
    ],
    inputs2: Annotated[
        list[Path] | None,
        typer.Option("--input2", help="R2 FASTQ file(s). Repeat to pair with --input."),
    ] = None,
    front: Annotated[
        str | None, typer.Option("--front", help="Cutadapt 5' adapter for R1.")
    ] = None,
    adapter: Annotated[
        str | None, typer.Option("--adapter", help="Cutadapt 3' adapter for R1.")
    ] = None,
    anywhere: Annotated[
        str | None, typer.Option("--anywhere", help="Cutadapt anywhere adapter for R1.")
    ] = None,
    front2: Annotated[
        str | None, typer.Option("--front2", help="Cutadapt 5' adapter for R2.")
    ] = None,
    adapter2: Annotated[
        str | None, typer.Option("--adapter2", help="Cutadapt 3' adapter for R2.")
    ] = None,
    anywhere2: Annotated[
        str | None, typer.Option("--anywhere2", help="Cutadapt anywhere adapter for R2.")
    ] = None,
    reads_per_file: Annotated[int, typer.Option("--reads-per-file", min=1)] = 1000,
    max_files: Annotated[int, typer.Option("--max-files", min=0)] = 16,
    max_mismatches: Annotated[int, typer.Option("--max-mismatches", min=0)] = 2,
    min_match_rate: Annotated[float, typer.Option("--min-match-rate", min=0.0, max=1.0)] = 0.8,
    mode: Annotated[str, typer.Option("--mode", help="warn or error.")] = "warn",
    output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
) -> None:
    """Check Cutadapt primer patterns against a deterministic FASTQ sample."""

    if inputs2 is not None and len(inputs2) != len(inputs):
        raise typer.BadParameter("--input2 must contain the same number of files as --input")
    files = [("R1", path) for path in inputs]
    files.extend(("R2", path) for path in (inputs2 or []))
    report = check_fastq_primers(
        files,
        cutadapt={
            "front": front,
            "adapter": adapter,
            "anywhere": anywhere,
            "front2": front2,
            "adapter2": adapter2,
            "anywhere2": anywhere2,
        },
        primer_check={
            "mode": mode,
            "reads_per_file": reads_per_file,
            "max_files": max_files,
            "max_mismatches": max_mismatches,
            "min_match_rate": min_match_rate,
        },
    )
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if output is not None:
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered + "\n", encoding="utf-8")
    typer.echo(rendered)
    if primer_check_fails(report, mode):
        raise typer.Exit(1)


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
