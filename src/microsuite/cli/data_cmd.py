from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.data.gut_to_soil import fetch_gut_to_soil
from microsuite.data.moving_pictures import copy_small_fixture, fetch_moving_pictures

app = typer.Typer(help="Fetch or materialize example datasets.", no_args_is_help=True)


@app.command("fetch")
def fetch(
    name: Annotated[str, typer.Argument(help="Dataset name.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output directory.")],
    full: Annotated[
        bool, typer.Option("--full", help="Download the fuller remote dataset.")
    ] = False,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing files.")] = False,
) -> None:
    if name == "gut-to-soil":
        fetch_gut_to_soil(output, force=force)
        return
    if name != "moving-pictures":
        raise typer.BadParameter("Available datasets: moving-pictures, gut-to-soil.")
    if full:
        fetch_moving_pictures(output, force=force)
    else:
        copy_small_fixture(output, force=force)
