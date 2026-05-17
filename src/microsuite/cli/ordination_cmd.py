from __future__ import annotations

from pathlib import Path
from typing import Annotated

import pandas as pd
import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.ordination.pcoa import pcoa

app = typer.Typer(help="Ordination commands.", no_args_is_help=True)


@app.command("pcoa")
def pcoa_cmd(
    distance_matrix: Annotated[Path, typer.Argument(help="Square TSV distance matrix.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output coordinates TSV.")],
    dimensions: Annotated[int, typer.Option("--dimensions", "-d", min=1)] = 3,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    dist = pd.read_csv(ensure_input(distance_matrix), sep="\t", index_col=0)
    result = pcoa(dist, dimensions=dimensions)
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)
