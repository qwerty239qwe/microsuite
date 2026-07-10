# src/microsuite/cli/refdb_cmd.py
from __future__ import annotations

from typing import Annotated

import typer

from microsuite.refdb.paths import VALID_BUILD_TARGETS
from microsuite.refdb.service import fetch_refdb
from microsuite.refdb.spec import RefDbSpec

app = typer.Typer(help="Fetch, build, and cache reference databases.", no_args_is_help=True)


@app.command("fetch")
def fetch(
    name: Annotated[str, typer.Argument(help="Reference DB name, e.g. homd, silva.")],
    version: Annotated[str, typer.Option("--version", help="DB version.")],
    provider: Annotated[str, typer.Option("--provider", help="Acquisition provider.")] = "biodbs",
    build: Annotated[
        str, typer.Option("--build", help="Build target: vsearch, blast, or qiime2.")
    ] = "vsearch",
    force: Annotated[bool, typer.Option("--force", help="Rebuild even if cached.")] = False,
) -> None:
    if build not in VALID_BUILD_TARGETS:
        raise typer.BadParameter(f"--build must be one of: {', '.join(VALID_BUILD_TARGETS)}")
    spec = RefDbSpec(name=name, version=version, provider=provider, build_targets=(build,))
    artifact = fetch_refdb(spec, build, force=force)
    typer.echo(str(artifact.path))
