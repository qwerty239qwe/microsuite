from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite._paths import ensure_input, prepare_output
from microsuite.io.biom import read_biom
from microsuite.io.h5ad import write_h5ad
from microsuite.io.qza import read_qza
from microsuite.io.tsv import read_tsv

app = typer.Typer(help="Import feature tables into AnnData .h5ad files.", no_args_is_help=True)


@app.command("tsv")
def import_tsv(
    table: Annotated[Path, typer.Argument(help="Feature table TSV.")],
    metadata: Annotated[Path, typer.Option("--metadata", "-m", help="Sample metadata TSV.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output .h5ad path.")],
    taxonomy: Annotated[Path | None, typer.Option("--taxonomy", "-t")] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_tsv(
        ensure_input(table),
        ensure_input(metadata),
        ensure_input(taxonomy) if taxonomy else None,
    )
    write_h5ad(adata, prepare_output(output, force=force))


@app.command("biom")
def import_biom(
    table: Annotated[Path, typer.Argument(help="BIOM feature table.")],
    metadata: Annotated[Path, typer.Option("--metadata", "-m", help="Sample metadata TSV.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output .h5ad path.")],
    taxonomy: Annotated[Path | None, typer.Option("--taxonomy", "-t")] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_biom(
        ensure_input(table),
        ensure_input(metadata),
        ensure_input(taxonomy) if taxonomy else None,
    )
    write_h5ad(adata, prepare_output(output, force=force))


@app.command("qza")
def import_qza(
    artifact: Annotated[Path, typer.Argument(help="QIIME2 .qza artifact.")],
    metadata: Annotated[Path, typer.Option("--metadata", "-m", help="Sample metadata TSV.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output .h5ad path.")],
    taxonomy_artifact: Annotated[
        Path | None, typer.Option("--taxonomy-artifact", help="Optional taxonomy .qza artifact.")
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    adata = read_qza(
        ensure_input(artifact),
        ensure_input(metadata),
        ensure_input(taxonomy_artifact) if taxonomy_artifact else None,
    )
    write_h5ad(adata, prepare_output(output, force=force))
