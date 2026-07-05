from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.cli._method_api import (
    evaluate,
    phylogeny,
    tax_barplot,
    tax_classify,
    tax_collapse,
    tax_train,
)


def register(app: typer.Typer) -> None:
    @app.command("tax_classify")
    def tax_classify_cmd(
        backend: Annotated[
            str,
            typer.Option(
                "--backend",
                "--method",
                help="Classification backend. --method is a deprecated alias.",
            ),
        ],
        rep_seqs: Annotated[
            Path,
            typer.Option("--rep-seqs", help="Representative sequences, usually a .qza artifact."),
        ],
        output: Annotated[Path, typer.Option("--output", "-o", help="Classification output.")],
        classifier: Annotated[
            str | None,
            typer.Option(
                "--classifier",
                help=(
                    "QIIME classifier, Kraken2/Bracken database, optional MetaPhlAn database, "
                    "or optional EMU database, as a path OR a cached reference as "
                    "'refdb:<name>@<version>[:<build>]'."
                ),
            ),
        ] = None,
        input_type: Annotated[
            str,
            typer.Option(
                "--input-type",
                help=(
                    "MetaPhlAn input type: fastq, fasta, bowtie2out, or sam. "
                    "EMU type: map-ont, map-pb, lr:hq, map-hifi, or sr; "
                    "default fastq maps to map-ont."
                ),
            ),
        ] = "fastq",
        level: Annotated[
            str,
            typer.Option("--level", help="Bracken taxonomy level: D, P, C, O, F, G, or S."),
        ] = "S",
        read_length: Annotated[
            int,
            typer.Option(
                "--read-length", min=1, help="Bracken read length used for abundance estimates."
            ),
        ] = 150,
        threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        tax_classify(
            backend=backend,
            rep_seqs=rep_seqs,
            classifier=classifier,
            input_type=input_type,
            level=level,
            read_length=read_length,
            output=output,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("tax_train")
    def tax_train_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        ref_seqs: Annotated[Path | None, typer.Option("--ref-seqs")] = None,
        ref_taxonomy: Annotated[Path | None, typer.Option("--ref-taxonomy")] = None,
        f_primer: Annotated[str | None, typer.Option("--f-primer")] = None,
        r_primer: Annotated[str | None, typer.Option("--r-primer")] = None,
        trunc_len: Annotated[int, typer.Option("--trunc-len", min=0)] = 0,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        threads: Annotated[str, typer.Option("--threads")] = "1",
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        tax_train(
            backend=backend,
            ref_seqs=ref_seqs,
            ref_taxonomy=ref_taxonomy,
            f_primer=f_primer,
            r_primer=r_primer,
            trunc_len=trunc_len,
            output=output,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("tax_barplot")
    def tax_barplot_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        table: Annotated[Path | None, typer.Option("--table")] = None,
        taxonomy: Annotated[Path | None, typer.Option("--taxonomy")] = None,
        metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        tax_barplot(
            backend=backend,
            table=table,
            taxonomy=taxonomy,
            metadata=metadata,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("tax_collapse")
    def tax_collapse_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        table: Annotated[Path | None, typer.Option("--table")] = None,
        taxonomy: Annotated[Path | None, typer.Option("--taxonomy")] = None,
        level: Annotated[int, typer.Option("--level", min=1)] = 6,
        output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        tax_collapse(
            backend=backend,
            table=table,
            taxonomy=taxonomy,
            level=level,
            output=output,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("phylogeny")
    def phylogeny_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Backend.")],
        rep_seqs: Annotated[Path | None, typer.Option("--rep-seqs")] = None,
        output_aligned: Annotated[Path | None, typer.Option("--output-aligned")] = None,
        output_masked: Annotated[Path | None, typer.Option("--output-masked")] = None,
        output_tree: Annotated[Path | None, typer.Option("--output-tree")] = None,
        output_rooted_tree: Annotated[Path | None, typer.Option("--output-rooted-tree")] = None,
        threads: Annotated[str, typer.Option("--threads")] = "1",
        force: Annotated[bool, typer.Option("--force")] = False,
        run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
        timeout: Annotated[float | None, typer.Option("--timeout")] = None,
    ) -> None:
        phylogeny(
            backend=backend,
            rep_seqs=rep_seqs,
            output_aligned=output_aligned,
            output_masked=output_masked,
            output_tree=output_tree,
            output_rooted_tree=output_rooted_tree,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("evaluate")
    def evaluate_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Evaluation backend.")],
        expected_taxa: Annotated[
            Path, typer.Option("--expected-taxa", help="Expected taxonomy artifact.")
        ],
        observed_taxa: Annotated[
            Path, typer.Option("--observed-taxa", help="Observed taxonomy artifact.")
        ],
        output: Annotated[Path, typer.Option("--output", "-o", help="Output visualization.")],
        feature_table: Annotated[
            Path | None, typer.Option("--feature-table", help="Optional feature table artifact.")
        ] = None,
        depth: Annotated[int, typer.Option("--depth", min=1)] = 7,
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        evaluate(
            backend=backend,
            expected_taxa=expected_taxa,
            observed_taxa=observed_taxa,
            feature_table=feature_table,
            output=output,
            depth=depth,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
