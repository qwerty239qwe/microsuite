from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.cli._method_api import (
    abundance,
    feature_filter,
    feature_summarize,
    normalize,
    rarefy,
    shared_taxa,
)


def register(app: typer.Typer) -> None:
    @app.command("normalize")
    def normalize_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Normalization backend.")],
        method: Annotated[
            str,
            typer.Option(
                "--method",
                help="relative, total-sum, clr, or prevalence-filter.",
            ),
        ],
        table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output .h5ad table.")],
        target_sum: Annotated[float, typer.Option("--target-sum")] = 1_000_000.0,
        pseudocount: Annotated[float, typer.Option("--pseudocount")] = 1.0,
        min_prevalence: Annotated[float, typer.Option("--min-prevalence")] = 0.1,
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    ) -> None:
        normalize(
            backend=backend,
            method=method,
            table=table,
            output=output,
            target_sum=target_sum,
            pseudocount=pseudocount,
            min_prevalence=min_prevalence,
            force=force,
        )

    @app.command("abundance")
    def abundance_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Abundance backend.")],
        table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
        level: Annotated[str, typer.Option("--level", help="Taxonomy level.")],
        relative: Annotated[bool, typer.Option("--relative/--counts")] = True,
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    ) -> None:
        abundance(
            backend=backend,
            table=table,
            output=output,
            level=level,
            relative=relative,
            force=force,
        )

    @app.command("shared_taxa")
    def shared_taxa_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Shared taxa backend.")],
        table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
        level: Annotated[str, typer.Option("--level", help="Taxonomy level.")],
        group: Annotated[str, typer.Option("--group", help="Sample metadata group column.")],
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    ) -> None:
        shared_taxa(
            backend=backend,
            table=table,
            output=output,
            level=level,
            group=group,
            force=force,
        )

    @app.command("rarefy")
    def rarefy_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Rarefaction backend.")],
        table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output .h5ad table.")],
        depth: Annotated[int, typer.Option("--depth", min=1, help="Rarefaction depth.")],
        seed: Annotated[int, typer.Option("--seed", help="Random seed.")] = 0,
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    ) -> None:
        rarefy(
            backend=backend,
            table=table,
            output=output,
            depth=depth,
            seed=seed,
            force=force,
        )

    @app.command("feature_summarize")
    def feature_summarize_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        mode: Annotated[str, typer.Option("--mode", help="summarize or tabulate-seqs")],
        table: Annotated[Path | None, typer.Option("--table")] = None,
        rep_seqs: Annotated[Path | None, typer.Option("--rep-seqs")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        feature_summarize(
            backend=backend,
            mode=mode,
            table=table,
            rep_seqs=rep_seqs,
            metadata=metadata,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("feature_filter")
    def feature_filter_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        table: Annotated[Path | None, typer.Option("--table")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        where: Annotated[str | None, typer.Option("--where")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        feature_filter(
            backend=backend,
            table=table,
            metadata=metadata,
            where=where,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
