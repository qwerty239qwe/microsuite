from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.methods.abundance import SUPPORTED_BACKENDS as ABUNDANCE_BACKENDS
from microsuite.methods.abundance import abundance
from microsuite.methods.cluster import SUPPORTED_BACKENDS as CLUSTER_BACKENDS
from microsuite.methods.cluster import cluster
from microsuite.methods.decontam import SUPPORTED_BACKENDS as DECONTAM_BACKENDS
from microsuite.methods.decontam import decontam
from microsuite.methods.denoise import SUPPORTED_BACKENDS as DENOISE_BACKENDS
from microsuite.methods.denoise import denoise
from microsuite.methods.diff_abundance import SUPPORTED_BACKENDS as DIFF_ABUNDANCE_BACKENDS
from microsuite.methods.diff_abundance import diff_abundance
from microsuite.methods.diversity_calc import SUPPORTED_METHODS as DIVERSITY_METHODS
from microsuite.methods.diversity_calc import diversity_calc
from microsuite.methods.evaluate import SUPPORTED_BACKENDS as EVALUATE_BACKENDS
from microsuite.methods.evaluate import evaluate
from microsuite.methods.functional_profile import SUPPORTED_BACKENDS as FUNCTIONAL_PROFILE_BACKENDS
from microsuite.methods.functional_profile import functional_profile
from microsuite.methods.network import SUPPORTED_BACKENDS as NETWORK_BACKENDS
from microsuite.methods.normalize import SUPPORTED_BACKENDS as NORMALIZE_BACKENDS
from microsuite.methods.normalize import normalize
from microsuite.methods.qc import SUPPORTED_BACKENDS as QC_BACKENDS
from microsuite.methods.qc import qc
from microsuite.methods.qc_filter import SUPPORTED_BACKENDS as QC_FILTER_BACKENDS
from microsuite.methods.qc_filter import qc_filter
from microsuite.methods.qiime2_wrappers import SUPPORTED_METHODS as QIIME2_WRAPPER_METHODS
from microsuite.methods.qiime2_wrappers import (
    demux,
    diff_viz,
    diversity_core,
    diversity_test,
    feature_filter,
    feature_summarize,
    metadata_tabulate,
    ordination_plot,
    phylogeny,
    qiime_import,
    rarefaction,
    tax_barplot,
    tax_collapse,
    tax_train,
)
from microsuite.methods.rarefy import SUPPORTED_BACKENDS as RAREFY_BACKENDS
from microsuite.methods.rarefy import rarefy
from microsuite.methods.report import SUPPORTED_BACKENDS as REPORT_BACKENDS
from microsuite.methods.report import report
from microsuite.methods.shared_taxa import SUPPORTED_BACKENDS as SHARED_TAXA_BACKENDS
from microsuite.methods.shared_taxa import shared_taxa
from microsuite.methods.tax_classify import SUPPORTED_METHODS, tax_classify
from microsuite.methods.trim import SUPPORTED_BACKENDS as TRIM_BACKENDS
from microsuite.methods.trim import trim

app = typer.Typer(help="Method-oriented microbiome operations.", no_args_is_help=True)


@app.command("methods")
def methods() -> None:
    typer.echo("qc")
    for backend in QC_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("qc_filter")
    for backend in QC_FILTER_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("trim")
    for backend in TRIM_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("denoise")
    for backend in DENOISE_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("cluster")
    for backend in CLUSTER_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("normalize")
    for backend in NORMALIZE_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("abundance")
    for backend in ABUNDANCE_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("shared_taxa")
    for backend in SHARED_TAXA_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("rarefy")
    for backend in RAREFY_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("tax_classify")
    for backend in SUPPORTED_METHODS:
        typer.echo(f"  - {backend}")
    typer.echo("diversity_calc")
    for backend in DIVERSITY_METHODS:
        typer.echo(f"  - {backend}")
    typer.echo("diff_abundance")
    for backend in DIFF_ABUNDANCE_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("decontam")
    for backend in DECONTAM_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("evaluate")
    for backend in EVALUATE_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("functional_profile")
    for backend in FUNCTIONAL_PROFILE_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("network")
    for backend in NETWORK_BACKENDS:
        typer.echo(f"  - {backend}")
    typer.echo("report")
    for backend in REPORT_BACKENDS:
        typer.echo(f"  - {backend}")
    for method, backends in QIIME2_WRAPPER_METHODS.items():
        typer.echo(method)
        for backend in backends:
            typer.echo(f"  - {backend}")


@app.command("qc")
def qc_cmd(
    backend: Annotated[str, typer.Option("--backend", help="QC backend.")],
    inputs: Annotated[
        list[Path] | None,
        typer.Option("--input", help="Input FASTQ file. Repeat for multiple files."),
    ] = None,
    input_dir: Annotated[
        Path | None, typer.Option("--input-dir", help="Input directory for MultiQC.")
    ] = None,
    demux: Annotated[
        Path | None, typer.Option("--demux", help="QIIME 2 demultiplexed reads artifact.")
    ] = None,
    output_dir: Annotated[
        Path | None, typer.Option("--output-dir", help="Output directory for report files.")
    ] = None,
    output: Annotated[
        Path | None, typer.Option("--output", "-o", help="Output QIIME 2 visualization.")
    ] = None,
    threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
    extract: Annotated[
        bool, typer.Option("--extract", help="Extract FastQC zip output after analysis.")
    ] = False,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
) -> None:
    qc(
        backend=backend,
        inputs=inputs,
        input_dir=input_dir,
        demux=demux,
        output_dir=output_dir,
        output=output,
        threads=threads,
        extract=extract,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


@app.command("metadata_tabulate")
def metadata_tabulate_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Backend.")],
    input_file: Annotated[Path | None, typer.Option("--input-file")] = None,
    output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
    force: Annotated[bool, typer.Option("--force")] = False,
    run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
    timeout: Annotated[float | None, typer.Option("--timeout")] = None,
) -> None:
    metadata_tabulate(
        backend=backend,
        input_file=input_file,
        output=output,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


@app.command("qiime_import")
def qiime_import_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Backend.")],
    input_path: Annotated[Path | None, typer.Option("--input-path")] = None,
    output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
    force: Annotated[bool, typer.Option("--force")] = False,
    run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
    timeout: Annotated[float | None, typer.Option("--timeout")] = None,
) -> None:
    qiime_import(
        backend=backend,
        input_path=input_path,
        output=output,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


@app.command("demux")
def demux_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Backend.")],
    seqs: Annotated[Path | None, typer.Option("--seqs")] = None,
    metadata: Annotated[Path | None, typer.Option("--metadata", "-m")] = None,
    barcode_column: Annotated[str | None, typer.Option("--barcode-column")] = None,
    output_demux: Annotated[Path | None, typer.Option("--output-demux")] = None,
    output_details: Annotated[Path | None, typer.Option("--output-details")] = None,
    force: Annotated[bool, typer.Option("--force")] = False,
    run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
    timeout: Annotated[float | None, typer.Option("--timeout")] = None,
) -> None:
    demux(
        backend=backend,
        seqs=seqs,
        metadata=metadata,
        barcode_column=barcode_column,
        output_demux=output_demux,
        output_details=output_details,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


@app.command("qc_filter")
def qc_filter_cmd(
    backend: Annotated[str, typer.Option("--backend", help="QC filtering backend.")],
    demux: Annotated[
        Path | None, typer.Option("--demux", help="QIIME 2 demultiplexed reads artifact.")
    ] = None,
    database: Annotated[
        Path | None, typer.Option("--database", help="QIIME 2 Bowtie2 index artifact.")
    ] = None,
    output: Annotated[
        Path | None,
        typer.Option(
            "--output",
            "-o",
            help="Output artifact for bowtie2-build or filtered demultiplexed reads.",
        ),
    ] = None,
    sequences: Annotated[
        Path | None,
        typer.Option("--sequences", help="Reference sequences for qiime2-bowtie2-build."),
    ] = None,
    query_sequences: Annotated[
        Path | None, typer.Option("--query-sequences", help="QIIME 2 query sequences artifact.")
    ] = None,
    reference_sequences: Annotated[
        Path | None,
        typer.Option("--reference-sequences", help="QIIME 2 reference sequences artifact."),
    ] = None,
    sequence_hits: Annotated[
        Path | None, typer.Option("--sequence-hits", help="Output sequences matching reference.")
    ] = None,
    sequence_misses: Annotated[
        Path | None,
        typer.Option("--sequence-misses", help="Output sequences not matching reference."),
    ] = None,
    method: Annotated[str, typer.Option("--method", help="Alignment method for exclude-seqs.")] = (
        "blast"
    ),
    perc_identity: Annotated[float, typer.Option("--perc-identity", min=0.0, max=1.0)] = 0.97,
    perc_query_aligned: Annotated[
        float, typer.Option("--perc-query-aligned", min=0.0, max=1.0)
    ] = 0.97,
    threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
    mode: Annotated[str, typer.Option("--mode", help="Bowtie2 alignment mode.")] = "local",
    sensitivity: Annotated[
        str, typer.Option("--sensitivity", help="Bowtie2 sensitivity preset.")
    ] = "sensitive",
    exclude: Annotated[bool, typer.Option("--exclude/--keep-matches")] = True,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
) -> None:
    qc_filter(
        backend=backend,
        demux=demux,
        database=database,
        output=output,
        sequences=sequences,
        query_sequences=query_sequences,
        reference_sequences=reference_sequences,
        sequence_hits=sequence_hits,
        sequence_misses=sequence_misses,
        method=method,
        perc_identity=perc_identity,
        perc_query_aligned=perc_query_aligned,
        threads=threads,
        mode=mode,
        sensitivity=sensitivity,
        exclude=exclude,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


@app.command("trim")
def trim_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Trim backend.")],
    read1: Annotated[Path, typer.Option("--read1", help="Forward or single-end FASTQ.")],
    output1: Annotated[Path, typer.Option("--output1", help="Output forward/single FASTQ.")],
    read2: Annotated[Path | None, typer.Option("--read2", help="Reverse FASTQ.")] = None,
    output2: Annotated[Path | None, typer.Option("--output2", help="Output reverse FASTQ.")] = None,
    unpaired1: Annotated[
        Path | None, typer.Option("--unpaired1", help="Trimmomatic unpaired R1 output.")
    ] = None,
    unpaired2: Annotated[
        Path | None, typer.Option("--unpaired2", help="Trimmomatic unpaired R2 output.")
    ] = None,
    html: Annotated[Path | None, typer.Option("--html", help="fastp HTML report.")] = None,
    json_report: Annotated[
        Path | None, typer.Option("--json-report", help="fastp or Cutadapt JSON report.")
    ] = None,
    adapter: Annotated[
        str | None, typer.Option("--adapter", "-a", help="Cutadapt 3' adapter for R1.")
    ] = None,
    front: Annotated[
        str | None, typer.Option("--front", "-g", help="Cutadapt 5' adapter for R1.")
    ] = None,
    anywhere: Annotated[
        str | None, typer.Option("--anywhere", "-b", help="Cutadapt 5'/3' adapter for R1.")
    ] = None,
    adapter2: Annotated[
        str | None, typer.Option("--adapter2", "-A", help="Cutadapt 3' adapter for R2.")
    ] = None,
    front2: Annotated[
        str | None, typer.Option("--front2", "-G", help="Cutadapt 5' adapter for R2.")
    ] = None,
    anywhere2: Annotated[
        str | None, typer.Option("--anywhere2", "-B", help="Cutadapt 5'/3' adapter for R2.")
    ] = None,
    quality_cutoff: Annotated[
        str | None,
        typer.Option("--quality-cutoff", "-q", help="Cutadapt quality cutoff, e.g. 20 or 10,20."),
    ] = None,
    minimum_length: Annotated[
        str | None, typer.Option("--minimum-length", "-m", help="Cutadapt minimum length.")
    ] = None,
    maximum_length: Annotated[
        str | None, typer.Option("--maximum-length", "-M", help="Cutadapt maximum length.")
    ] = None,
    max_n: Annotated[
        str | None, typer.Option("--max-n", help="Discard reads with more than this many Ns.")
    ] = None,
    discard_untrimmed: Annotated[
        bool, typer.Option("--discard-untrimmed", help="Discard reads without adapter matches.")
    ] = False,
    trimmomatic_steps: Annotated[
        list[str] | None,
        typer.Option(
            "--trimmomatic-step",
            help="Trimmomatic trimming step. Repeat for multiple steps.",
        ),
    ] = None,
    basename: Annotated[
        str | None, typer.Option("--basename", help="Trim Galore output basename.")
    ] = None,
    trim_galore_version: Annotated[
        str,
        typer.Option(
            "--trim-galore-version",
            help="Trim Galore compatibility mode: auto, legacy, or v2.",
        ),
    ] = "auto",
    threads: Annotated[str, typer.Option("--threads", help="Thread count or 'auto'.")] = "1",
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
) -> None:
    trim(
        backend=backend,
        read1=read1,
        read2=read2,
        output1=output1,
        output2=output2,
        unpaired1=unpaired1,
        unpaired2=unpaired2,
        html=html,
        json_report=json_report,
        adapter=adapter,
        front=front,
        anywhere=anywhere,
        adapter2=adapter2,
        front2=front2,
        anywhere2=anywhere2,
        quality_cutoff=quality_cutoff,
        minimum_length=minimum_length,
        maximum_length=maximum_length,
        max_n=max_n,
        discard_untrimmed=discard_untrimmed,
        trimmomatic_steps=trimmomatic_steps,
        basename=basename,
        trim_galore_version=trim_galore_version,
        threads=threads,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


@app.command("denoise")
def denoise_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Denoising backend.")],
    demux: Annotated[Path, typer.Option("--demux", help="Demultiplexed reads artifact.")],
    output_table: Annotated[
        Path, typer.Option("--output-table", help="Output feature table artifact.")
    ],
    output_rep_seqs: Annotated[
        Path, typer.Option("--output-rep-seqs", help="Output representative sequences artifact.")
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
    paired: Annotated[bool, typer.Option("--paired", help="Use paired-end DADA2 mode.")] = False,
    trim_left: Annotated[int, typer.Option("--trim-left", min=0)] = 0,
    trunc_len: Annotated[int, typer.Option("--trunc-len", min=0)] = 0,
    trim_left_f: Annotated[int, typer.Option("--trim-left-f", min=0)] = 0,
    trunc_len_f: Annotated[int, typer.Option("--trunc-len-f", min=0)] = 0,
    trim_left_r: Annotated[int, typer.Option("--trim-left-r", min=0)] = 0,
    trunc_len_r: Annotated[int, typer.Option("--trunc-len-r", min=0)] = 0,
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
        paired=paired,
        trim_left=trim_left,
        trunc_len=trunc_len,
        trim_left_f=trim_left_f,
        trunc_len_f=trunc_len_f,
        trim_left_r=trim_left_r,
        trunc_len_r=trunc_len_r,
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
        float, typer.Option("--identity", min=0.0, max=1.0, help="Clustering identity threshold.")
    ] = 0.97,
    table: Annotated[
        Path | None,
        typer.Option("--table", help="Input feature table artifact. Required for qiime2-vsearch."),
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
        Path, typer.Option("--rep-seqs", help="Representative sequences, usually a .qza artifact.")
    ],
    output: Annotated[Path, typer.Option("--output", "-o", help="Classification output.")],
    classifier: Annotated[
        Path | None,
        typer.Option(
            "--classifier",
            help=(
                "QIIME classifier, Kraken2/Bracken database, optional MetaPhlAn database, "
                "or optional EMU database."
            ),
        ),
    ] = None,
    input_type: Annotated[
        str,
        typer.Option(
            "--input-type",
            help=(
                "MetaPhlAn input type: fastq, fasta, bowtie2out, or sam. "
                "EMU type: map-ont, map-pb, lr:hq, map-hifi, or sr; default fastq maps to map-ont."
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


@app.command("diff_viz")
def diff_viz_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Backend.")],
    data: Annotated[Path | None, typer.Option("--data")] = None,
    output: Annotated[Path | None, typer.Option("--output", "-o")] = None,
    force: Annotated[bool, typer.Option("--force")] = False,
    run_dir: Annotated[Path | None, typer.Option("--run-dir")] = None,
    timeout: Annotated[float | None, typer.Option("--timeout")] = None,
) -> None:
    diff_viz(
        backend=backend,
        data=data,
        output=output,
        force=force,
        run_dir=run_dir,
        timeout=timeout,
    )


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


@app.command("diff_abundance")
def diff_abundance_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Differential abundance backend.")],
    table: Annotated[Path, typer.Option("--table", help="Input .h5ad table.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output TSV.")],
    group: Annotated[str, typer.Option("--group", help="Sample metadata group column.")],
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
    )


@app.command("decontam")
def decontam_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Contaminant detection backend.")],
    table: Annotated[Path, typer.Option("--table", help="QIIME 2 feature table artifact.")],
    metadata: Annotated[Path, typer.Option("--metadata", help="Sample metadata TSV.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output decontam score artifact.")],
    method: Annotated[
        str, typer.Option("--method", help="prevalence, frequency, or combined.")
    ] = "prevalence",
    prev_control_column: Annotated[
        str | None, typer.Option("--prev-control-column", help="Negative-control metadata column.")
    ] = None,
    prev_control_indicator: Annotated[
        str | None,
        typer.Option("--prev-control-indicator", help="Negative-control metadata value."),
    ] = None,
    freq_concentration_column: Annotated[
        str | None,
        typer.Option("--freq-concentration-column", help="Sample concentration metadata column."),
    ] = None,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
    run_dir: Annotated[
        Path | None, typer.Option("--run-dir", help="Write runtime logs here.")
    ] = None,
    timeout: Annotated[
        float | None, typer.Option("--timeout", help="Command timeout in seconds.")
    ] = None,
) -> None:
    decontam(
        backend=backend,
        table=table,
        metadata=metadata,
        output=output,
        method=method,
        prev_control_column=prev_control_column,
        prev_control_indicator=prev_control_indicator,
        freq_concentration_column=freq_concentration_column,
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
