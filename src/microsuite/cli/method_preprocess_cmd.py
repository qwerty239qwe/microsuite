from __future__ import annotations

from pathlib import Path
from typing import Annotated

import typer

from microsuite.cli._method_api import (
    decontam,
    qc,
    qc_filter,
    trim,
)


def register(app: typer.Typer) -> None:
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
            Path | None,
            typer.Option("--sequence-hits", help="Output sequences matching reference."),
        ] = None,
        sequence_misses: Annotated[
            Path | None,
            typer.Option("--sequence-misses", help="Output sequences not matching reference."),
        ] = None,
        method: Annotated[
            str, typer.Option("--method", help="Alignment method for exclude-seqs.")
        ] = ("blast"),
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
        output2: Annotated[
            Path | None, typer.Option("--output2", help="Output reverse FASTQ.")
        ] = None,
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
            typer.Option(
                "--quality-cutoff", "-q", help="Cutadapt quality cutoff, e.g. 20 or 10,20."
            ),
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

    @app.command("decontam")
    def decontam_cmd(
        backend: Annotated[str, typer.Option("--backend", help="Contaminant detection backend.")],
        table: Annotated[Path, typer.Option("--table", help="QIIME 2 feature table artifact.")],
        metadata: Annotated[Path, typer.Option("--metadata", help="Sample metadata TSV.")],
        output: Annotated[
            Path, typer.Option("--output", "-o", help="Output decontam score artifact.")
        ],
        method: Annotated[
            str, typer.Option("--method", help="prevalence, frequency, or combined.")
        ] = "prevalence",
        prev_control_column: Annotated[
            str | None,
            typer.Option("--prev-control-column", help="Negative-control metadata column."),
        ] = None,
        prev_control_indicator: Annotated[
            str | None,
            typer.Option("--prev-control-indicator", help="Negative-control metadata value."),
        ] = None,
        freq_concentration_column: Annotated[
            str | None,
            typer.Option(
                "--freq-concentration-column", help="Sample concentration metadata column."
            ),
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
