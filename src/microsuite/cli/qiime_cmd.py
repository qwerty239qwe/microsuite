from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Annotated

import typer

from microsuite._paths import ensure_input
from microsuite.qiime2.artifact import extract_data_payload, inspect_artifact

app = typer.Typer(help="QIIME 2 artifact utilities.", no_args_is_help=True)


@app.command("inspect")
def inspect(
    artifact: Annotated[Path, typer.Argument(help="QIIME 2 .qza or .qzv artifact.")],
    json_output: Annotated[bool, typer.Option("--json", help="Emit JSON.")] = False,
) -> None:
    info = inspect_artifact(ensure_input(artifact))
    if json_output:
        typer.echo(json.dumps(asdict(info), indent=2))
        return
    typer.echo(f"path: {info.path}")
    typer.echo(f"uuid: {info.uuid}")
    typer.echo(f"type: {info.artifact_type}")
    typer.echo(f"format: {info.format}")
    typer.echo(f"framework: {info.framework_version}")
    typer.echo(f"archive: {info.archive_version}")
    typer.echo("data:")
    for data_file in info.data_files:
        typer.echo(f"  - {data_file}")


@app.command("extract")
def extract(
    artifact: Annotated[Path, typer.Argument(help="QIIME 2 .qza or .qzv artifact.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output directory.")],
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing files.")] = False,
) -> None:
    written = extract_data_payload(ensure_input(artifact), output, force=force)
    for path in written:
        typer.echo(path)
