from __future__ import annotations

from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods._qiime import ensure_inputs, prepare_outputs, require_qiime, run_qiime
from microsuite.runtime.runner import resolve_threads

SUPPORTED_BACKENDS = (
    "qiime2-quality-filter-q-score",
    "qiime2-bowtie2-build",
    "qiime2-filter-reads",
    "qiime2-exclude-seqs",
)


def qc_filter(
    *,
    backend: str,
    demux: Path | None = None,
    database: Path | None = None,
    output: Path | None = None,
    sequences: Path | None = None,
    query_sequences: Path | None = None,
    reference_sequences: Path | None = None,
    sequence_hits: Path | None = None,
    sequence_misses: Path | None = None,
    method: str = "blast",
    perc_identity: float = 0.97,
    perc_query_aligned: float = 0.97,
    threads: int | str = 1,
    mode: str = "local",
    sensitivity: str = "sensitive",
    exclude: bool = True,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    resolved_threads = resolve_threads(threads)
    if backend == "qiime2-bowtie2-build":
        qc_filter_qiime2_bowtie2_build(
            sequences=sequences,
            output=output,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "qiime2-quality-filter-q-score":
        qc_filter_qiime2_quality_filter_q_score(
            demux=demux,
            output=output,
            sequence_hits=sequence_hits,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "qiime2-filter-reads":
        qc_filter_qiime2_filter_reads(
            demux=demux,
            database=database,
            output=output,
            threads=resolved_threads,
            mode=mode,
            sensitivity=sensitivity,
            exclude=exclude,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "qiime2-exclude-seqs":
        qc_filter_qiime2_exclude_seqs(
            query_sequences=query_sequences,
            reference_sequences=reference_sequences,
            sequence_hits=sequence_hits,
            sequence_misses=sequence_misses,
            method=method,
            perc_identity=perc_identity,
            perc_query_aligned=perc_query_aligned,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    raise MicrobiomeSuiteError(
        f"Unsupported QC filter backend '{backend}'. Choose one of: {', '.join(SUPPORTED_BACKENDS)}"
    )


def qc_filter_qiime2_quality_filter_q_score(
    *,
    demux: Path | None,
    output: Path | None,
    sequence_hits: Path | None,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if demux is None:
        raise MicrobiomeSuiteError(
            "--demux is required for --backend qiime2-quality-filter-q-score."
        )
    if output is None:
        raise MicrobiomeSuiteError(
            "--output is required for --backend qiime2-quality-filter-q-score."
        )
    if sequence_hits is None:
        raise MicrobiomeSuiteError(
            "--sequence-hits is required for --backend qiime2-quality-filter-q-score."
        )
    qiime = require_qiime("QIIME 2 quality-filter q-score")
    ensure_inputs(demux)
    prepare_outputs(output, sequence_hits, force=force)
    command = [
        qiime,
        "quality-filter",
        "q-score",
        "--i-demux",
        str(demux),
        "--o-filtered-sequences",
        str(output),
        "--o-filter-stats",
        str(sequence_hits),
    ]
    run_qiime(
        command,
        "QIIME 2 quality-filter q-score failed.",
        run_dir=run_dir,
        timeout=timeout,
        task="qc_filter",
        backend="qiime2-quality-filter-q-score",
    )


def qc_filter_qiime2_bowtie2_build(
    *,
    sequences: Path | None,
    output: Path | None,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if sequences is None:
        raise MicrobiomeSuiteError("--sequences is required for --backend qiime2-bowtie2-build.")
    if output is None:
        raise MicrobiomeSuiteError("--output is required for --backend qiime2-bowtie2-build.")
    qiime = require_qiime("QIIME 2 quality-control bowtie2-build")
    ensure_inputs(sequences)
    prepare_outputs(output, force=force)

    command = [
        qiime,
        "quality-control",
        "bowtie2-build",
        "--i-sequences",
        str(sequences),
        "--p-n-threads",
        str(threads),
        "--o-database",
        str(output),
    ]
    run_qiime(
        command,
        "QIIME 2 quality-control bowtie2-build failed.",
        run_dir=run_dir,
        timeout=timeout,
        task="qc_filter",
        backend="qiime2-bowtie2-build",
    )


def qc_filter_qiime2_filter_reads(
    *,
    demux: Path | None,
    database: Path | None,
    output: Path | None,
    threads: int,
    mode: str,
    sensitivity: str,
    exclude: bool,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if demux is None:
        raise MicrobiomeSuiteError("--demux is required for --backend qiime2-filter-reads.")
    if database is None:
        raise MicrobiomeSuiteError("--database is required for --backend qiime2-filter-reads.")
    if output is None:
        raise MicrobiomeSuiteError("--output is required for --backend qiime2-filter-reads.")
    qiime = require_qiime("QIIME 2 quality-control filter-reads")
    ensure_inputs(demux, database)
    prepare_outputs(output, force=force)

    command = [
        qiime,
        "quality-control",
        "filter-reads",
        "--i-demultiplexed-sequences",
        str(demux),
        "--i-database",
        str(database),
        "--p-n-threads",
        str(threads),
        "--p-mode",
        mode,
        "--p-sensitivity",
        sensitivity,
        "--p-exclude-seqs" if exclude else "--p-no-exclude-seqs",
        "--o-filtered-sequences",
        str(output),
    ]
    run_qiime(
        command,
        "QIIME 2 quality-control filter-reads failed.",
        run_dir=run_dir,
        timeout=timeout,
        task="qc_filter",
        backend="qiime2-filter-reads",
    )


def qc_filter_qiime2_exclude_seqs(
    *,
    query_sequences: Path | None,
    reference_sequences: Path | None,
    sequence_hits: Path | None,
    sequence_misses: Path | None,
    method: str,
    perc_identity: float,
    perc_query_aligned: float,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if query_sequences is None:
        raise MicrobiomeSuiteError(
            "--query-sequences is required for --backend qiime2-exclude-seqs."
        )
    if reference_sequences is None:
        raise MicrobiomeSuiteError(
            "--reference-sequences is required for --backend qiime2-exclude-seqs."
        )
    if sequence_hits is None:
        raise MicrobiomeSuiteError("--sequence-hits is required for --backend qiime2-exclude-seqs.")
    if sequence_misses is None:
        raise MicrobiomeSuiteError(
            "--sequence-misses is required for --backend qiime2-exclude-seqs."
        )
    qiime = require_qiime("QIIME 2 quality-control exclude-seqs")
    ensure_inputs(query_sequences, reference_sequences)
    prepare_outputs(sequence_hits, sequence_misses, force=force)

    command = [
        qiime,
        "quality-control",
        "exclude-seqs",
        "--i-query-sequences",
        str(query_sequences),
        "--i-reference-sequences",
        str(reference_sequences),
        "--p-method",
        method,
        "--p-perc-identity",
        str(perc_identity),
        "--p-perc-query-aligned",
        str(perc_query_aligned),
    ]
    if method == "vsearch":
        command.extend(["--p-threads", str(threads)])
    command.extend(
        [
            "--o-sequence-hits",
            str(sequence_hits),
            "--o-sequence-misses",
            str(sequence_misses),
        ]
    )
    run_qiime(
        command,
        "QIIME 2 quality-control exclude-seqs failed.",
        run_dir=run_dir,
        timeout=timeout,
        task="qc_filter",
        backend="qiime2-exclude-seqs",
    )
