from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.examples.moving_pictures import run_example

app = typer.Typer(help="Run complete example pipelines.", no_args_is_help=True)


@app.command("moving-pictures")
def moving_pictures(
    output: Annotated[Path, typer.Option("--out", "-o", help="Output run directory.")],
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
) -> None:
    run_example(output, force=force)
