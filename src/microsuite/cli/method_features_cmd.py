from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.cli._method_api import (
    assemble,
    bin_contigs,
    cluster,
    denoise,
)


def register(app: typer.Typer) -> None:
    @app.command("denoise")
    def denoise_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Denoising backend.")],
        demux: Annotated[Path, typer.Option("--demux", help="Demultiplexed reads artifact.")],
        output_table: Annotated[
            Path, typer.Option("--output-table", help="Output feature table artifact.")
        ],
        output_rep_seqs: Annotated[
            Path,
            typer.Option("--output-rep-seqs", help="Output representative sequences artifact."),
        ],
        output_stats: Annotated[
            Path, typer.Option("--output-stats", help="Output denoising stats artifact.")
        ],
        output_base_transition_stats: Annotated[
            Path | None,
            typer.Option(
                "--output-base-transition-stats",
                help="Optional DADA2 base transition stats.",
            ),
        ] = None,
        output_base_transition_plot: Annotated[
            Path | None,
            typer.Option(
                "--output-base-transition-plot",
                help=(
                    "Optional DADA2 base transition visualization (.qzv). "
                    "Requires --output-base-transition-stats."
                ),
            ),
        ] = None,
        mode: Annotated[
            str | None,
            typer.Option("--mode", help="DADA2 mode: single, paired, ccs, or pyro."),
        ] = None,
        paired: Annotated[
            bool,
            typer.Option("--paired", help="Deprecated alias for --mode paired."),
        ] = False,
        trim_left: Annotated[int, typer.Option("--trim-left", min=0)] = 0,
        trunc_len: Annotated[int, typer.Option("--trunc-len", min=0)] = 0,
        trim_left_f: Annotated[int, typer.Option("--trim-left-f", min=0)] = 0,
        trunc_len_f: Annotated[int, typer.Option("--trunc-len-f", min=0)] = 0,
        trim_left_r: Annotated[int, typer.Option("--trim-left-r", min=0)] = 0,
        trunc_len_r: Annotated[int, typer.Option("--trunc-len-r", min=0)] = 0,
        max_ee: Annotated[
            float | None,
            typer.Option("--max-ee", help="DADA2 single/CCS/R maximum expected errors."),
        ] = None,
        max_ee_f: Annotated[
            float | None,
            typer.Option("--max-ee-f", help="DADA2 paired forward maximum expected errors."),
        ] = None,
        max_ee_r: Annotated[
            float | None,
            typer.Option("--max-ee-r", help="DADA2 paired reverse maximum expected errors."),
        ] = None,
        trunc_q: Annotated[
            int | None,
            typer.Option("--trunc-q", min=0, help="DADA2 quality-score truncation threshold."),
        ] = None,
        max_n: Annotated[
            int | None, typer.Option("--max-n", min=0, help="R/DADA2 maximum ambiguous bases.")
        ] = None,
        rm_phix: Annotated[
            bool | None, typer.Option("--rm-phix/--no-rm-phix", help="R/DADA2 PhiX removal toggle.")
        ] = None,
        pooling_method: Annotated[
            str | None,
            typer.Option("--pooling-method", help="DADA2 pooling method: independent or pseudo."),
        ] = None,
        chimera_method: Annotated[
            str | None,
            typer.Option("--chimera-method", help="DADA2 chimera method: consensus or none."),
        ] = None,
        min_fold_parent_over_abundance: Annotated[
            float | None,
            typer.Option(
                "--min-fold-parent-over-abundance",
                help="DADA2 chimera parent abundance fold threshold.",
            ),
        ] = None,
        allow_one_off: Annotated[
            bool | None,
            typer.Option(
                "--allow-one-off/--no-allow-one-off",
                help="DADA2 one-off chimera detection toggle.",
            ),
        ] = None,
        n_reads_learn: Annotated[
            int | None,
            typer.Option(
                "--n-reads-learn", min=1, help="Reads used to train the DADA2 error model."
            ),
        ] = None,
        hashed_feature_ids: Annotated[
            bool | None,
            typer.Option(
                "--hashed-feature-ids/--no-hashed-feature-ids",
                help="QIIME DADA2 hashed feature IDs.",
            ),
        ] = None,
        retain_all_samples: Annotated[
            bool | None,
            typer.Option(
                "--retain-all-samples/--drop-empty-samples",
                help="QIIME DADA2 sample retention toggle.",
            ),
        ] = None,
        min_overlap: Annotated[
            int | None,
            typer.Option("--min-overlap", min=4, help="Paired DADA2 minimum merge overlap."),
        ] = None,
        max_merge_mismatch: Annotated[
            int | None,
            typer.Option(
                "--max-merge-mismatch", min=0, help="Paired DADA2 maximum merge mismatches."
            ),
        ] = None,
        trim_overhang: Annotated[
            bool | None,
            typer.Option(
                "--trim-overhang/--no-trim-overhang",
                help="Paired DADA2 merge overhang trimming toggle.",
            ),
        ] = None,
        ccs_front: Annotated[
            str | None, typer.Option("--ccs-front", help="QIIME DADA2 CCS front primer sequence.")
        ] = None,
        ccs_adapter: Annotated[
            str | None, typer.Option("--ccs-adapter", help="QIIME DADA2 CCS adapter sequence.")
        ] = None,
        ccs_max_mismatch: Annotated[
            int | None,
            typer.Option(
                "--ccs-max-mismatch",
                min=0,
                help="QIIME DADA2 CCS primer/adapter mismatches.",
            ),
        ] = None,
        ccs_indels: Annotated[
            bool | None,
            typer.Option(
                "--ccs-indels/--no-ccs-indels",
                help="QIIME DADA2 CCS primer/adapter indel matching.",
            ),
        ] = None,
        ccs_min_len: Annotated[
            int | None,
            typer.Option("--ccs-min-len", min=0, help="QIIME DADA2 CCS minimum read length."),
        ] = None,
        ccs_max_len: Annotated[
            int | None,
            typer.Option("--ccs-max-len", min=0, help="QIIME DADA2 CCS maximum read length."),
        ] = None,
        threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        denoise(
            backend=backend,
            demux=demux,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            output_stats=output_stats,
            output_base_transition_stats=output_base_transition_stats,
            output_base_transition_plot=output_base_transition_plot,
            mode=mode,
            paired=paired,
            trim_left=trim_left,
            trunc_len=trunc_len,
            trim_left_f=trim_left_f,
            trunc_len_f=trunc_len_f,
            trim_left_r=trim_left_r,
            trunc_len_r=trunc_len_r,
            max_ee=max_ee,
            max_ee_f=max_ee_f,
            max_ee_r=max_ee_r,
            trunc_q=trunc_q,
            max_n=max_n,
            rm_phix=rm_phix,
            pooling_method=pooling_method,
            chimera_method=chimera_method,
            min_fold_parent_over_abundance=min_fold_parent_over_abundance,
            allow_one_off=allow_one_off,
            n_reads_learn=n_reads_learn,
            hashed_feature_ids=hashed_feature_ids,
            retain_all_samples=retain_all_samples,
            min_overlap=min_overlap,
            max_merge_mismatch=max_merge_mismatch,
            trim_overhang=trim_overhang,
            ccs_front=ccs_front,
            ccs_adapter=ccs_adapter,
            ccs_max_mismatch=ccs_max_mismatch,
            ccs_indels=ccs_indels,
            ccs_min_len=ccs_min_len,
            ccs_max_len=ccs_max_len,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("cluster")
    def cluster_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Clustering backend.")],
        rep_seqs: Annotated[
            Path, typer.Option("--rep-seqs", help="Input representative sequences artifact.")
        ],
        output_table: Annotated[
            Path, typer.Option("--output-table", help="Output clustered table artifact.")
        ],
        output_rep_seqs: Annotated[
            Path, typer.Option("--output-rep-seqs", help="Output clustered sequences artifact.")
        ],
        identity: Annotated[
            float,
            typer.Option("--identity", min=0.0, max=1.0, help="Clustering identity threshold."),
        ] = 0.97,
        table: Annotated[
            Path | None,
            typer.Option(
                "--table", help="Input feature table artifact. Required for qiime2-vsearch."
            ),
        ] = None,
        output_uc: Annotated[
            Path | None,
            typer.Option("--output-uc", help="Optional USEARCH/VSEARCH .uc cluster mapping."),
        ] = None,
        sample_delimiter: Annotated[
            str,
            typer.Option("--sample-delimiter", help="Delimiter for sample IDs in sequence labels."),
        ] = "_",
        sample_field: Annotated[
            int,
            typer.Option("--sample-field", min=0, help="Zero-based sample ID field in labels."),
        ] = 0,
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        cluster(
            backend=backend,
            rep_seqs=rep_seqs,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            table=table,
            identity=identity,
            output_uc=output_uc,
            sample_delimiter=sample_delimiter,
            sample_field=sample_field,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("assemble")
    def assemble_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Assembly backend.")],
        output_dir: Annotated[
            Path | None,
            typer.Option("--output-dir", help="Output directory."),
        ] = None,
        read1: Annotated[
            Path | None, typer.Option("--read1", help="Forward or single-end FASTQ.")
        ] = None,
        read2: Annotated[Path | None, typer.Option("--read2", help="Reverse FASTQ.")] = None,
        reads: Annotated[
            Path | None,
            typer.Option("--reads", help="Single/interleaved reads file. IDBA-UD expects FASTA."),
        ] = None,
        output_contigs: Annotated[
            Path | None,
            typer.Option("--output-contigs", help="MOSHPIT output contigs artifact."),
        ] = None,
        presets: Annotated[
            str,
            typer.Option("--presets", help="MOSHPIT MEGAHIT preset."),
        ] = "meta-sensitive",
        min_contig: Annotated[
            int,
            typer.Option("--min-contig", min=1, help="MOSHPIT minimum contig length."),
        ] = 500,
        parallel_config: Annotated[
            Path | None,
            typer.Option("--parallel-config", help="MOSHPIT parsl parallel config TOML."),
        ] = None,
        verbose: Annotated[
            bool, typer.Option("--verbose", help="Pass --verbose to MOSHPIT.")
        ] = False,
        threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        assemble(
            backend=backend,
            read1=read1,
            read2=read2,
            reads=reads,
            output_dir=output_dir,
            output_contigs=output_contigs,
            presets=presets,
            min_contig=min_contig,
            parallel_config=parallel_config,
            verbose=verbose,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )

    @app.command("bin")
    def bin_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Binning backend.")],
        contigs: Annotated[Path, typer.Option("--contigs", help="Input contigs FASTA.")],
        output_dir: Annotated[
            Path | None,
            typer.Option("--output-dir", help="Output directory."),
        ] = None,
        depth: Annotated[
            Path | None, typer.Option("--depth", help="MetaBAT2 depth matrix.")
        ] = None,
        abundance: Annotated[
            Path | None, typer.Option("--abundance", help="MaxBin2 abundance table.")
        ] = None,
        coverage: Annotated[
            Path | None, typer.Option("--coverage", help="CONCOCT coverage table.")
        ] = None,
        alignment_maps: Annotated[
            Path | None,
            typer.Option("--alignment-maps", help="MOSHPIT read-to-contig alignment maps."),
        ] = None,
        output_mags: Annotated[
            Path | None, typer.Option("--output-mags", help="MOSHPIT output MAGs artifact.")
        ] = None,
        output_contig_map: Annotated[
            Path | None,
            typer.Option("--output-contig-map", help="MOSHPIT output contig-map artifact."),
        ] = None,
        output_unbinned_contigs: Annotated[
            Path | None,
            typer.Option(
                "--output-unbinned-contigs",
                help="MOSHPIT output unbinned-contigs artifact.",
            ),
        ] = None,
        prefix: Annotated[str, typer.Option("--prefix", help="Output bin prefix.")] = "bin",
        seed: Annotated[
            int, typer.Option("--seed", min=0, help="MOSHPIT random seed.")
        ] = 100,
        parallel_config: Annotated[
            Path | None,
            typer.Option("--parallel-config", help="MOSHPIT parsl parallel config TOML."),
        ] = None,
        verbose: Annotated[
            bool, typer.Option("--verbose", help="Pass --verbose to MOSHPIT.")
        ] = False,
        threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
        force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
        run_dir: Annotated[
            Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
        ] = None,
        timeout: Annotated[
            float | None, typer.Option("--timeout", help="Command timeout in seconds.")
        ] = None,
    ) -> None:
        bin_contigs(
            backend=backend,
            contigs=contigs,
            depth=depth,
            abundance=abundance,
            coverage=coverage,
            alignment_maps=alignment_maps,
            output_mags=output_mags,
            output_contig_map=output_contig_map,
            output_unbinned_contigs=output_unbinned_contigs,
            output_dir=output_dir,
            prefix=prefix,
            seed=seed,
            parallel_config=parallel_config,
            verbose=verbose,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
