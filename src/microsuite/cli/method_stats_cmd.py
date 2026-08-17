from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.cli._method_api import (
    diff_abundance,
    functional_profile,
    report,
)


def register(app: typer.Typer) -> None:
    @app.command("diff_abundance")
    def diff_abundance_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Differential abundance backend.")],
        table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
        output: Annotated[
            Path,
            typer.Option("--output", "-o", help="Output TSV, or output directory for MaAsLin 3."),
        ],
        group: Annotated[
            str | None,
            typer.Option("--group", help="Sample metadata group column."),
        ] = None,
        metadata: Annotated[
            Path | None, typer.Option("--metadata", "-m", help="QIIME 2 sample metadata TSV.")
        ] = None,
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
        runtime: Annotated[
            str,
            typer.Option("--runtime", help="R backend runtime: 'local' Rscript or 'docker'."),
        ] = "local",
        image: Annotated[
            str | None,
            typer.Option("--image", help="Override the per-backend r-diffab container image."),
        ] = None,
        formula: Annotated[
            str | None,
            typer.Option(
                "--formula",
                help="Complete MaAsLin 3 lme4 formula, including random effects if needed.",
            ),
        ] = None,
        fix_formula: Annotated[
            str | None,
            typer.Option("--fix-formula", help="MaAsLin 3 fixed-effects formula RHS."),
        ] = None,
        rand_formula: Annotated[
            str | None,
            typer.Option(
                "--rand-formula", help="MaAsLin 3 random-effects term, e.g. '(1|subject)'."
            ),
        ] = None,
        normalization: Annotated[
            str,
            typer.Option("--normalization", help="MaAsLin 3: TSS, CLR, or NONE."),
        ] = "TSS",
        transform: Annotated[
            str,
            typer.Option("--transform", help="MaAsLin 3: LOG, PLOG, or NONE."),
        ] = "LOG",
        min_prevalence: Annotated[
            float,
            typer.Option("--min-prevalence", help="MaAsLin 3 feature prevalence cutoff."),
        ] = 0.0,
        min_abundance: Annotated[
            float,
            typer.Option("--min-abundance", help="MaAsLin 3 feature abundance cutoff."),
        ] = 0.0,
    ) -> None:
        diff_abundance(
            backend=backend,
            table=table,
            group=group,
            output=output,
            metadata=metadata,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
            runtime=runtime,
            image=image,
            formula=formula,
            fix_formula=fix_formula,
            rand_formula=rand_formula,
            normalization=normalization,
            transform=transform,
            min_prevalence=min_prevalence,
            min_abundance=min_abundance,
        )

    @app.command("functional_profile")
    def functional_profile_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Functional profiling backend.")],
        output_dir: Annotated[Path, typer.Option("--output-dir", help="Output directory.")],
        table: Annotated[Path | None, typer.Option("--table", help="Feature/OTU table.")] = None,
        rep_seqs: Annotated[
            Path | None,
            typer.Option("--rep-seqs", help="Representative sequences FASTA."),
        ] = None,
        reads: Annotated[
            Path | None,
            typer.Option("--reads", help="Metagenomic reads for HUMAnN."),
        ] = None,
        database: Annotated[
            Path | None,
            typer.Option("--database", help="Reference database directory."),
        ] = None,
        protein_database: Annotated[
            Path | None,
            typer.Option("--protein-database", help="HUMAnN protein database directory."),
        ] = None,
        threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
        database_mode: Annotated[
            str,
            typer.Option("--database-mode", help="Tax4Fun2 database mode: Ref99NR or Ref100NR."),
        ] = "Ref99NR",
        min_identity: Annotated[
            float,
            typer.Option("--min-identity", min=0.0, max=1.0, help="Tax4Fun2 identity threshold."),
        ] = 0.97,
        normalize_pathways: Annotated[
            bool,
            typer.Option("--normalize-pathways", help="Normalize Tax4Fun2 pathway assignments."),
        ] = False,
        force: Annotated[
            bool, typer.Option("--force", help="Allow non-empty output directory.")
        ] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        functional_profile(
            backend=backend,
            table=table,
            rep_seqs=rep_seqs,
            reads=reads,
            database=database,
            protein_database=protein_database,
            output_dir=output_dir,
            threads=threads,
            database_mode=database_mode,
            min_identity=min_identity,
            normalize_pathways=normalize_pathways,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("report")
    def report_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Report backend.")],
        run_dir: Annotated[Path, typer.Option("--run-dir", help="Input run directory.")],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output HTML report.")],
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    ) -> None:
        report(backend=backend, run_dir=run_dir, output=output, force=force)
