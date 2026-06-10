from __future__ import annotations

import shutil
from importlib.resources import files
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("qiime2-dada2", "qiime2-deblur", "dada2-r")
DADA2_MODES = ("single", "paired", "ccs", "pyro")
POOLING_METHODS = ("independent", "pseudo")
CHIMERA_METHODS = ("consensus", "none")
DADA2_R_SCRIPT = "dada2_denoise.R"


def denoise(
    *,
    backend: str,
    demux: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    output_base_transition_stats: Path | None = None,
    output_base_transition_plot: Path | None = None,
    mode: str | None = None,
    paired: bool = False,
    trim_left: int = 0,
    trunc_len: int = 0,
    trim_left_f: int = 0,
    trunc_len_f: int = 0,
    trim_left_r: int = 0,
    trunc_len_r: int = 0,
    max_ee: float | None = None,
    max_ee_f: float | None = None,
    max_ee_r: float | None = None,
    trunc_q: int | None = None,
    max_n: int | None = None,
    rm_phix: bool | None = None,
    pooling_method: str | None = None,
    chimera_method: str | None = None,
    min_fold_parent_over_abundance: float | None = None,
    allow_one_off: bool | None = None,
    n_reads_learn: int | None = None,
    hashed_feature_ids: bool | None = None,
    retain_all_samples: bool | None = None,
    min_overlap: int | None = None,
    max_merge_mismatch: int | None = None,
    trim_overhang: bool | None = None,
    ccs_front: str | None = None,
    ccs_adapter: str | None = None,
    ccs_max_mismatch: int | None = None,
    ccs_indels: bool | None = None,
    ccs_min_len: int | None = None,
    ccs_max_len: int | None = None,
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
            output_base_transition_plot=output_base_transition_plot,
            mode=_resolve_dada2_mode(mode=mode, paired=paired),
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
        dada2_r_mode = _resolve_dada2_mode(mode=mode, paired=paired)
        if dada2_r_mode not in ("single", "paired"):
            raise MicrobiomeSuiteError("R/DADA2 denoising supports only single or paired mode.")
        denoise_dada2_r(
            input_dir=demux,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            output_stats=output_stats,
            paired=dada2_r_mode == "paired",
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
            min_overlap=min_overlap,
            max_merge_mismatch=max_merge_mismatch,
            trim_overhang=trim_overhang,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
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
    output_base_transition_plot: Path | None,
    mode: str,
    trim_left: int,
    trunc_len: int,
    trim_left_f: int,
    trunc_len_f: int,
    trim_left_r: int,
    trunc_len_r: int,
    max_ee: float | None,
    max_ee_f: float | None,
    max_ee_r: float | None,
    trunc_q: int | None,
    pooling_method: str | None,
    chimera_method: str | None,
    min_fold_parent_over_abundance: float | None,
    allow_one_off: bool | None,
    n_reads_learn: int | None,
    hashed_feature_ids: bool | None,
    retain_all_samples: bool | None,
    min_overlap: int | None,
    max_merge_mismatch: int | None,
    trim_overhang: bool | None,
    ccs_front: str | None,
    ccs_adapter: str | None,
    ccs_max_mismatch: int | None,
    ccs_indels: bool | None,
    ccs_min_len: int | None,
    ccs_max_len: int | None,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    _validate_dada2_common(pooling_method=pooling_method, chimera_method=chimera_method)
    if output_base_transition_plot is not None and output_base_transition_stats is None:
        raise MicrobiomeSuiteError(
            "--output-base-transition-plot requires --output-base-transition-stats."
        )
    if mode != "paired" and any(
        value is not None for value in (min_overlap, max_merge_mismatch, trim_overhang)
    ):
        raise MicrobiomeSuiteError(
            "--min-overlap, --max-merge-mismatch, and --trim-overhang only apply "
            "to DADA2 paired mode."
        )
    if mode != "ccs" and any(
        value is not None
        for value in (
            ccs_front,
            ccs_adapter,
            ccs_max_mismatch,
            ccs_indels,
            ccs_min_len,
            ccs_max_len,
        )
    ):
        raise MicrobiomeSuiteError("--ccs-* options only apply to DADA2 ccs mode.")
    if mode == "ccs" and ccs_front is None:
        raise MicrobiomeSuiteError("--ccs-front is required for DADA2 ccs mode.")

    qiime = _require_qiime("QIIME 2 DADA2 denoising")
    ensure_input(demux)
    _prepare_outputs(output_table, output_rep_seqs, output_stats, force=force)
    if output_base_transition_stats is not None:
        _prepare_outputs(output_base_transition_stats, force=force)
    if output_base_transition_plot is not None:
        _prepare_outputs(output_base_transition_plot, force=force)

    command = [qiime, "dada2", f"denoise-{mode}"]
    command.extend(["--i-demultiplexed-seqs", str(demux)])
    if mode == "paired":
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
        _append_value(command, "--p-max-ee-f", max_ee_f)
        _append_value(command, "--p-max-ee-r", max_ee_r)
        _append_value(command, "--p-min-overlap", min_overlap)
        _append_value(command, "--p-max-merge-mismatch", max_merge_mismatch)
        _append_bool(command, "--p-trim-overhang", trim_overhang)
    elif mode == "ccs":
        command.extend(["--p-trunc-len", str(trunc_len), "--p-trim-left", str(trim_left)])
        _append_value(command, "--p-front", ccs_front)
        _append_value(command, "--p-adapter", ccs_adapter)
        _append_value(command, "--p-max-mismatch", ccs_max_mismatch)
        _append_bool(command, "--p-indels", ccs_indels)
        _append_value(command, "--p-min-len", ccs_min_len)
        _append_value(command, "--p-max-len", ccs_max_len)
        _append_value(command, "--p-max-ee", max_ee)
    else:
        command.extend(
            [
                "--p-trim-left",
                str(trim_left),
                "--p-trunc-len",
                str(trunc_len),
            ]
        )
        _append_value(command, "--p-max-ee", max_ee)
    _append_value(command, "--p-trunc-q", trunc_q)
    _append_value(command, "--p-pooling-method", pooling_method)
    _append_value(command, "--p-chimera-method", chimera_method)
    _append_value(
        command,
        "--p-min-fold-parent-over-abundance",
        min_fold_parent_over_abundance,
    )
    _append_bool(command, "--p-allow-one-off", allow_one_off)
    _append_value(command, "--p-n-reads-learn", n_reads_learn)
    _append_bool(command, "--p-hashed-feature-ids", hashed_feature_ids)
    _append_bool(command, "--p-retain-all-samples", retain_all_samples)
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
    params = _dada2_log_params(
        mode=mode,
        max_ee=max_ee,
        max_ee_f=max_ee_f,
        max_ee_r=max_ee_r,
        trunc_q=trunc_q,
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
    )
    _run(
        command,
        "QIIME 2 DADA2 denoising failed.",
        run_dir=run_dir,
        timeout=timeout,
        backend="qiime2-dada2",
        params=params,
    )
    if output_base_transition_plot is not None:
        plot_command = [
            qiime,
            "dada2",
            "plot-base-transitions",
            "--i-base-transition-stats",
            str(output_base_transition_stats),
            "--o-visualization",
            str(output_base_transition_plot),
        ]
        _run(
            plot_command,
            "QIIME 2 DADA2 base-transition plotting failed.",
            run_dir=run_dir,
            timeout=timeout,
            backend="qiime2-dada2",
            params={**params, "base_transition_plot": str(output_base_transition_plot)},
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
    max_ee: float | None,
    max_ee_f: float | None,
    max_ee_r: float | None,
    trunc_q: int | None,
    max_n: int | None,
    rm_phix: bool | None,
    pooling_method: str | None,
    chimera_method: str | None,
    min_fold_parent_over_abundance: float | None,
    allow_one_off: bool | None,
    n_reads_learn: int | None,
    min_overlap: int | None,
    max_merge_mismatch: int | None,
    trim_overhang: bool | None,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    _validate_dada2_common(pooling_method=pooling_method, chimera_method=chimera_method)
    if not paired and any(
        value is not None for value in (min_overlap, max_merge_mismatch, trim_overhang)
    ):
        raise MicrobiomeSuiteError(
            "--min-overlap, --max-merge-mismatch, and --trim-overhang only apply "
            "to paired DADA2 R mode."
        )
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
        _append_value(command, "--max-ee-f", max_ee_f)
        _append_value(command, "--max-ee-r", max_ee_r)
        _append_value(command, "--min-overlap", min_overlap)
        _append_value(command, "--max-merge-mismatch", max_merge_mismatch)
        _append_bool(command, "--trim-overhang", trim_overhang)
    else:
        command.extend(["--trim-left", str(trim_left), "--trunc-len", str(trunc_len)])
        _append_value(command, "--max-ee", max_ee)
    _append_value(command, "--trunc-q", trunc_q)
    _append_value(command, "--max-n", max_n)
    _append_bool(command, "--rm-phix", rm_phix)
    _append_value(command, "--pooling-method", pooling_method)
    _append_value(command, "--chimera-method", chimera_method)
    _append_value(
        command,
        "--min-fold-parent-over-abundance",
        min_fold_parent_over_abundance,
    )
    _append_bool(command, "--allow-one-off", allow_one_off)
    _append_value(command, "--n-reads-learn", n_reads_learn)
    _run(
        command,
        "R/DADA2 denoising failed.",
        run_dir=run_dir,
        timeout=timeout,
        backend="dada2-r",
        params=_dada2_log_params(
            mode="paired" if paired else "single",
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
            min_overlap=min_overlap,
            max_merge_mismatch=max_merge_mismatch,
            trim_overhang=trim_overhang,
        ),
    )


def _resolve_dada2_mode(*, mode: str | None, paired: bool) -> str:
    if mode is None:
        return "paired" if paired else "single"
    normalized = mode.lower()
    if normalized not in DADA2_MODES:
        raise MicrobiomeSuiteError(
            f"Unsupported DADA2 mode '{mode}'. Choose one of: {', '.join(DADA2_MODES)}"
        )
    if paired and normalized != "paired":
        raise MicrobiomeSuiteError("--paired is a deprecated alias for --mode paired.")
    return normalized


def _validate_dada2_common(
    *, pooling_method: str | None, chimera_method: str | None
) -> None:
    if pooling_method is not None and pooling_method not in POOLING_METHODS:
        raise MicrobiomeSuiteError(
            f"Unsupported DADA2 pooling method '{pooling_method}'. "
            f"Choose one of: {', '.join(POOLING_METHODS)}"
        )
    if chimera_method is not None and chimera_method not in CHIMERA_METHODS:
        raise MicrobiomeSuiteError(
            f"Unsupported DADA2 chimera method '{chimera_method}'. "
            f"Choose one of: {', '.join(CHIMERA_METHODS)}"
        )


def _append_value(command: list[str], flag: str, value: object | None) -> None:
    if value is not None:
        command.extend([flag, str(value)])


def _append_bool(command: list[str], flag: str, value: bool | None) -> None:
    if value is not None:
        if value:
            command.append(flag)
        elif flag.startswith("--p-"):
            command.append(flag.replace("--p-", "--p-no-", 1))
        else:
            command.append(flag.replace("--", "--no-", 1))


def _dada2_log_params(**params: object | None) -> dict[str, object]:
    return {key: value for key, value in params.items() if value is not None}


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
    params: dict[str, object] | None = None,
) -> None:
    run_command(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="denoise", backend=backend, params=params or {}),
    )
