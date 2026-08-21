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
        subclass: Annotated[
            str | None,
            typer.Option("--subclass", help="LEfSe subclass/block metadata column."),
        ] = None,
        reference: Annotated[
            str | None,
            typer.Option(
                "--reference",
                help=(
                    "LEfSe reference class, or MaAsLin 3 'column,level' pairs separated "
                    "by ';'. Defaults to the first sorted level."
                ),
            ),
        ] = None,
        seed: Annotated[int, typer.Option("--seed", help="LEfSe random seed.")] = 1234,
        kruskal_threshold: Annotated[
            float, typer.Option("--kruskal-threshold", help="LEfSe Kruskal-Wallis p-value cutoff.")
        ] = 0.05,
        wilcoxon_threshold: Annotated[
            float, typer.Option("--wilcoxon-threshold", help="LEfSe Wilcoxon p-value cutoff.")
        ] = 0.05,
        lda_threshold: Annotated[
            float, typer.Option("--lda-threshold", help="LEfSe minimum absolute LDA score.")
        ] = 2.0,
        p_adjust_method: Annotated[
            str, typer.Option("--p-adjust-method", help="LEfSe p-value adjustment method.")
        ] = "none",
        trim_names: Annotated[
            bool, typer.Option("--trim-names", help="Trim LEfSe feature names to terminal labels.")
        ] = False,
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
            subclass=subclass,
            reference=reference,
            seed=seed,
            kruskal_threshold=kruskal_threshold,
            wilcoxon_threshold=wilcoxon_threshold,
            lda_threshold=lda_threshold,
            p_adjust_method=p_adjust_method,
            trim_names=trim_names,
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
        runtime: Annotated[
            str,
            typer.Option(
                "--runtime",
                help=(
                    "PICRUSt2/Tax4Fun2 execution: local executable(s) or the selected Docker image."
                ),
            ),
        ] = "local",
        image: Annotated[
            str | None,
            typer.Option(
                "--image",
                help="Override the PICRUSt2 or Tax4Fun2 functional-profile image.",
            ),
        ] = None,
        engine: Annotated[
            str,
            typer.Option("--engine", help="Container engine for --runtime docker."),
        ] = "docker",
        picrust2_database: Annotated[
            str,
            typer.Option(
                "--picrust2-database",
                help="PICRUSt2 database: SC, oldIMG, or custom (case-insensitive).",
            ),
        ] = "SC",
        picrust2_ref_dir1: Annotated[
            Path | None,
            typer.Option(
                "--picrust2-ref-dir1", help="Custom PICRUSt2 EPA-ng reference directory 1."
            ),
        ] = None,
        picrust2_ref_dir2: Annotated[
            Path | None,
            typer.Option(
                "--picrust2-ref-dir2", help="Custom PICRUSt2 EPA-ng reference directory 2."
            ),
        ] = None,
        picrust2_custom_trait_tables_ref1: Annotated[
            list[Path] | None,
            typer.Option(
                "--picrust2-custom-trait-tables-ref1",
                "--picrust2-custom-trait-table-ref1",
                "--picrust2-custom-trait-tables",
                help="Custom trait table for reference 1; repeat or comma-delimit.",
            ),
        ] = None,
        picrust2_custom_trait_tables_ref2: Annotated[
            list[Path] | None,
            typer.Option(
                "--picrust2-custom-trait-tables-ref2",
                "--picrust2-custom-trait-table-ref2",
                help="Custom trait table for reference 2; repeat or comma-delimit.",
            ),
        ] = None,
        picrust2_marker_gene_table_ref1: Annotated[
            Path | None,
            typer.Option(
                "--picrust2-marker-gene-table-ref1",
                "--picrust2-marker-gene-table",
                help="Custom marker-gene table for reference 1.",
            ),
        ] = None,
        picrust2_marker_gene_table_ref2: Annotated[
            Path | None,
            typer.Option(
                "--picrust2-marker-gene-table-ref2",
                help="Custom marker-gene table for reference 2.",
            ),
        ] = None,
        picrust2_pathway_map: Annotated[
            Path | None,
            typer.Option("--picrust2-pathway-map", help="PICRUSt2 pathway map file."),
        ] = None,
        picrust2_reaction_func: Annotated[
            str | None,
            typer.Option(
                "--picrust2-reaction-func", help="PICRUSt2 symbolic trait or reaction-map path."
            ),
        ] = None,
        picrust2_regroup_map: Annotated[
            Path | None,
            typer.Option("--picrust2-regroup-map", help="PICRUSt2 regroup map file."),
        ] = None,
        picrust2_no_regroup: Annotated[
            bool,
            typer.Option(
                "--picrust2-no-regroup",
                help="Skip PICRUSt2 regrouping and do not apply a regroup map.",
            ),
        ] = False,
        picrust2_no_pathways: Annotated[
            bool,
            typer.Option("--picrust2-no-pathways", help="Do not infer PICRUSt2 pathways."),
        ] = False,
        picrust2_coverage: Annotated[
            bool,
            typer.Option(
                "--picrust2-coverage", help="Request experimental PICRUSt2 pathway coverage."
            ),
        ] = False,
        picrust2_max_nsti: Annotated[
            float,
            typer.Option(
                "--picrust2-max-nsti", min=0.0, help="Maximum PICRUSt2 NSTI (default: 2.0)."
            ),
        ] = 2.0,
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
            runtime=runtime,
            image=image,
            engine=engine,
            picrust2_database=picrust2_database,
            picrust2_ref_dir1=picrust2_ref_dir1,
            picrust2_ref_dir2=picrust2_ref_dir2,
            picrust2_custom_trait_tables_ref1=picrust2_custom_trait_tables_ref1,
            picrust2_custom_trait_tables_ref2=picrust2_custom_trait_tables_ref2,
            picrust2_marker_gene_table_ref1=picrust2_marker_gene_table_ref1,
            picrust2_marker_gene_table_ref2=picrust2_marker_gene_table_ref2,
            picrust2_pathway_map=picrust2_pathway_map,
            picrust2_reaction_func=picrust2_reaction_func,
            picrust2_regroup_map=picrust2_regroup_map,
            picrust2_no_regroup=picrust2_no_regroup,
            picrust2_no_pathways=picrust2_no_pathways,
            picrust2_coverage=picrust2_coverage,
            picrust2_max_nsti=picrust2_max_nsti,
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
