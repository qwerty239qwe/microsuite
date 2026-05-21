from __future__ import annotations

import shutil
from importlib.resources import files
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("qiime2-dada2", "qiime2-deblur", "dada2-r")
PLANNED_BACKENDS: tuple[str, ...] = ()
DADA2_R_SCRIPT = "dada2_denoise.R"


def denoise(
    *,
    backend: str,
    demux: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    output_base_transition_stats: Path | None = None,
    paired: bool = False,
    trim_left: int = 0,
    trunc_len: int = 0,
    trim_left_f: int = 0,
    trunc_len_f: int = 0,
    trim_left_r: int = 0,
    trunc_len_r: int = 0,
    threads: int | str = 1,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    resolved_threads = resolve_threads(threads)
    if backend == "qiime2-dada2":
        denoise_qiime2_dada2(
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
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "qiime2-deblur":
        denoise_qiime2_deblur(
            demux=demux,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            output_stats=output_stats,
            trim_left=trim_left,
            trunc_len=trunc_len,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "dada2-r":
        denoise_dada2_r(
            input_dir=demux,
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
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend in PLANNED_BACKENDS:
        raise MicrobiomeSuiteError(
            f"Denoise backend '{backend}' is registered but not implemented yet. "
            "Use --backend qiime2-dada2 or --backend qiime2-deblur for now."
        )
    raise MicrobiomeSuiteError(
        f"Unsupported denoise backend '{backend}'. Choose one of: {', '.join(SUPPORTED_BACKENDS)}"
    )


def denoise_qiime2_dada2(
    *,
    demux: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    output_base_transition_stats: Path | None,
    paired: bool,
    trim_left: int,
    trunc_len: int,
    trim_left_f: int,
    trunc_len_f: int,
    trim_left_r: int,
    trunc_len_r: int,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    qiime = _require_qiime("QIIME 2 DADA2 denoising")
    ensure_input(demux)
    _prepare_outputs(output_table, output_rep_seqs, output_stats, force=force)
    if output_base_transition_stats is not None:
        _prepare_outputs(output_base_transition_stats, force=force)

    command = [qiime, "dada2", "denoise-paired" if paired else "denoise-single"]
    command.extend(["--i-demultiplexed-seqs", str(demux)])
    if paired:
        command.extend(
            [
                "--p-trim-left-f",
                str(trim_left_f),
                "--p-trunc-len-f",
                str(trunc_len_f),
                "--p-trim-left-r",
                str(trim_left_r),
                "--p-trunc-len-r",
                str(trunc_len_r),
            ]
        )
    else:
        command.extend(
            [
                "--p-trim-left",
                str(trim_left),
                "--p-trunc-len",
                str(trunc_len),
            ]
        )
    command.extend(
        [
            "--o-table",
            str(output_table),
            "--o-representative-sequences",
            str(output_rep_seqs),
            "--o-denoising-stats",
            str(output_stats),
            "--p-n-threads",
            str(threads),
        ]
    )
    if output_base_transition_stats is not None:
        command.extend(["--o-base-transition-stats", str(output_base_transition_stats)])
    _run(
        command,
        "QIIME 2 DADA2 denoising failed.",
        run_dir=run_dir,
        timeout=timeout,
        backend="qiime2-dada2",
    )


def denoise_qiime2_deblur(
    *,
    demux: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    trim_left: int,
    trunc_len: int,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if trunc_len < 1:
        raise MicrobiomeSuiteError(
            "--trunc-len must be greater than zero for --backend qiime2-deblur."
        )
    qiime = _require_qiime("QIIME 2 Deblur denoising")
    ensure_input(demux)
    _prepare_outputs(output_table, output_rep_seqs, output_stats, force=force)

    command = [
        qiime,
        "deblur",
        "denoise-16S",
        "--i-demultiplexed-seqs",
        str(demux),
        "--p-trim-length",
        str(trunc_len),
        "--p-left-trim-len",
        str(trim_left),
        "--p-jobs-to-start",
        str(threads),
        "--o-table",
        str(output_table),
        "--o-representative-sequences",
        str(output_rep_seqs),
        "--o-stats",
        str(output_stats),
    ]
    _run(
        command,
        "QIIME 2 Deblur denoising failed.",
        run_dir=run_dir,
        timeout=timeout,
        backend="qiime2-deblur",
    )


def denoise_dada2_r(
    *,
    input_dir: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    paired: bool,
    trim_left: int,
    trunc_len: int,
    trim_left_f: int,
    trunc_len_f: int,
    trim_left_r: int,
    trunc_len_r: int,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise MicrobiomeSuiteError(
            "R/DADA2 denoising requires the external 'Rscript' command. "
            "Install R with the dada2 package and rerun this command."
        )
    if not input_dir.exists() or not input_dir.is_dir():
        raise MicrobiomeSuiteError(f"Input directory does not exist: {input_dir}")
    _prepare_outputs(output_table, output_rep_seqs, output_stats, force=force)

    command = [
        rscript,
        str(files("microsuite.resources").joinpath(DADA2_R_SCRIPT)),
        "--input-dir",
        str(input_dir),
        "--output-table",
        str(output_table),
        "--output-rep-seqs",
        str(output_rep_seqs),
        "--output-stats",
        str(output_stats),
        "--threads",
        str(threads),
    ]
    if paired:
        command.append("--paired")
        command.extend(
            [
                "--trim-left-f",
                str(trim_left_f),
                "--trunc-len-f",
                str(trunc_len_f),
                "--trim-left-r",
                str(trim_left_r),
                "--trunc-len-r",
                str(trunc_len_r),
            ]
        )
    else:
        command.extend(["--trim-left", str(trim_left), "--trunc-len", str(trunc_len)])
    _run(command, "R/DADA2 denoising failed.", run_dir=run_dir, timeout=timeout, backend="dada2-r")


def _require_qiime(task: str) -> str:
    qiime = shutil.which("qiime")
    if qiime is None:
        raise MicrobiomeSuiteError(
            f"{task} requires the external 'qiime' command. "
            "Activate a QIIME 2 environment and rerun this command."
        )
    return qiime


def _prepare_outputs(*outputs: Path, force: bool) -> None:
    for output in outputs:
        prepare_output(output, force=force)


def _run(
    command: list[str],
    failure_message: str,
    *,
    run_dir: Path | None,
    timeout: float | None,
    backend: str,
) -> None:
    run_command(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="denoise", backend=backend),
    )
