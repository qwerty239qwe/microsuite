from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.cli._method_api import (
    diversity_calc,
    diversity_core,
    diversity_test,
    ordination_plot,
    rarefaction,
)


def register(app: typer.Typer) -> None:
    @app.command("diversity_calc")
    def diversity_calc_cmd(
        backend: Annotated[
            str,
            typer.Option(
                "--backend",
                "--method",
                help="Diversity backend. --method is a deprecated alias.",
            ),
        ],
        metric: Annotated[
            str,
            typer.Option(
                "--metric",
                help="Metric, e.g. shannon, bray-curtis, faith-pd, weighted-unifrac.",
            ),
        ],
        table: Annotated[Path, typer.Option("--table", help="QIIME 2 feature table artifact.")],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output QIIME 2 artifact.")],
        phylogeny: Annotated[
            Path | None, typer.Option("--phylogeny", help="Required for phylogenetic metrics.")
        ] = None,
        threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        diversity_calc(
            backend=backend,
            metric=metric,
            table=table,
            phylogeny=phylogeny,
            output=output,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("diversity_core")
    def diversity_core_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        table: Annotated[Path | None, typer.Option("--table")] = None,
        phylogeny_path: Annotated[Path | None, typer.Option("--phylogeny")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        sampling_depth: Annotated[int, typer.Option("--sampling-depth", min=1)] = 1103,
        output_dir: Annotated[Path | None, typer.Option("--output-dir")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        diversity_core(
            backend=backend,
            table=table,
            phylogeny_path=phylogeny_path,
            metadata=metadata,
            sampling_depth=sampling_depth,
            output_dir=output_dir,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("diversity_test")
    def diversity_test_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        alpha_diversity: Annotated[Path | None, typer.Option("--alpha-diversity")] = None,
        distance_matrix: Annotated[Path | None, typer.Option("--distance-matrix")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        metadata_column: Annotated[str | None, typer.Option("--metadata-column")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        method: Annotated[str, typer.Option("--method")] = "permanova",
        pairwise: Annotated[bool, typer.Option("--pairwise")] = False,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        diversity_test(
            backend=backend,
            alpha_diversity=alpha_diversity,
            distance_matrix=distance_matrix,
            metadata=metadata,
            metadata_column=metadata_column,
            output=output,
            method=method,
            pairwise=pairwise,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("rarefaction")
    def rarefaction_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        table: Annotated[Path | None, typer.Option("--table")] = None,
        phylogeny_path: Annotated[Path | None, typer.Option("--phylogeny")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        max_depth: Annotated[int, typer.Option("--max-depth", min=1)] = 4000,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        rarefaction(
            backend=backend,
            table=table,
            phylogeny_path=phylogeny_path,
            metadata=metadata,
            max_depth=max_depth,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("ordination_plot")
    def ordination_plot_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        pcoa: Annotated[Path | None, typer.Option("--pcoa")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        ordination_plot(
            backend=backend,
            pcoa=pcoa,
            metadata=metadata,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
