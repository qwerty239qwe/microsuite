from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.methods.abundance import SUPPORTED_BACKENDS as ABUNDANCE_BACKENDS
from microsuite.methods.abundance import abundance
from microsuite.methods.cluster import SUPPORTED_BACKENDS as CLUSTER_BACKENDS
from microsuite.methods.cluster import cluster
from microsuite.methods.denoise import SUPPORTED_BACKENDS as DENOISE_BACKENDS
from microsuite.methods.denoise import denoise
from microsuite.methods.diff_abundance import SUPPORTED_BACKENDS as DIFF_ABUNDANCE_BACKENDS
from microsuite.methods.diff_abundance import diff_abundance
from microsuite.methods.diversity_calc import SUPPORTED_METHODS as DIVERSITY_METHODS
from microsuite.methods.diversity_calc import diversity_calc
from microsuite.methods.normalize import SUPPORTED_BACKENDS as NORMALIZE_BACKENDS
from microsuite.methods.normalize import normalize
from microsuite.methods.rarefy import SUPPORTED_BACKENDS as RAREFY_BACKENDS
from microsuite.methods.rarefy import rarefy
from microsuite.methods.report import SUPPORTED_BACKENDS as REPORT_BACKENDS
from microsuite.methods.report import report
from microsuite.methods.shared_taxa import SUPPORTED_BACKENDS as SHARED_TAXA_BACKENDS
from microsuite.methods.shared_taxa import shared_taxa
from microsuite.methods.tax_classify import SUPPORTED_METHODS, tax_classify

app = typer.Typer(help="Method-oriented microbiome operations.", no_args_is_help=True)


@app.command("methods")
def methods() -> None:
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
    typer.echo("report")
    for backend in REPORT_BACKENDS:
        typer.echo(f"  - {backend}")


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
    paired: Annotated[bool, typer.Option("--paired", help="Use paired-end DADA2 mode.")] = False,
    trim_left: Annotated[int, typer.Option("--trim-left", min=0)] = 0,
    trunc_len: Annotated[int, typer.Option("--trunc-len", min=0)] = 0,
    trim_left_f: Annotated[int, typer.Option("--trim-left-f", min=0)] = 0,
    trunc_len_f: Annotated[int, typer.Option("--trunc-len-f", min=0)] = 0,
    trim_left_r: Annotated[int, typer.Option("--trim-left-r", min=0)] = 0,
    trunc_len_r: Annotated[int, typer.Option("--trunc-len-r", min=0)] = 0,
    threads: Annotated[int, typer.Option("--threads", min=1)] = 1,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
) -> None:
    denoise(
        backend=backend,
        demux=demux,
        output_table=output_table,
        output_rep_seqs=output_rep_seqs,
        output_stats=output_stats,
        paired=paired,
        trim_left=trim_left,
        trunc_len=trunc_len,
        trim_left_f=trim_left_f,
        trunc_len_f=trunc_len_f,
        trim_left_r=trim_left_r,
        trunc_len_r=trunc_len_r,
        threads=threads,
        force=force,
    )


@app.command("cluster")
def cluster_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Clustering backend.")],
    table: Annotated[Path, typer.Option("--table", help="Input feature table artifact.")],
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
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing outputs.")] = False,
) -> None:
    cluster(
        backend=backend,
        table=table,
        rep_seqs=rep_seqs,
        output_table=output_table,
        output_rep_seqs=output_rep_seqs,
        identity=identity,
        force=force,
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
        Path | None, typer.Option("--classifier", help="Pretrained classifier for qiime2.")
    ] = None,
    threads: Annotated[int, typer.Option("--threads", min=1)] = 1,
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    tax_classify(
        backend=backend,
        rep_seqs=rep_seqs,
        classifier=classifier,
        output=output,
        threads=threads,
        force=force,
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
) -> None:
    diversity_calc(
        backend=backend,
        metric=metric,
        table=table,
        phylogeny=phylogeny,
        output=output,
        threads=threads,
        force=force,
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
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    diff_abundance(
        backend=backend,
        table=table,
        group=group,
        output=output,
        force=force,
    )


@app.command("report")
def report_cmd(
    backend: Annotated[str, typer.Option("--backend", help="Report backend.")],
    run_dir: Annotated[Path, typer.Option("--run-dir", help="Input run directory.")],
    output: Annotated[Path, typer.Option("--output", "-o", help="Output HTML report.")],
    force: Annotated[bool, typer.Option("--force", help="Overwrite existing output.")] = False,
) -> None:
    report(backend=backend, run_dir=run_dir, output=output, force=force)
