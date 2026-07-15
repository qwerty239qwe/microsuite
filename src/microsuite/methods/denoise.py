from __future__ import annotations

import json
import re
import shutil
import warnings
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib.resources import files
from pathlib import Path

from microsuite import __version__ as _MICROSUITE_VERSION
from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.metadata import Artifact, ProvenanceFile, StageRecord, stage_execution
from microsuite.methods import dada2_manifest, dada2_qc
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command
from microsuite.runtime.validation import validate_outputs

SUPPORTED_BACKENDS = ("qiime2-dada2", "qiime2-deblur", "dada2-r")
DADA2_MODES = ("single", "paired", "ccs", "pyro")
POOLING_METHODS = ("independent", "pseudo")
CHIMERA_METHODS = ("consensus", "none")
DADA2_R_SCRIPT = "dada2_denoise.R"

_FASTQ_EXTS = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
_READ_PATTERNS = [
    re.compile(r"^(?P<sample>.+?)[._-]R(?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)[._-]read(?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)[._-](?P<read>[12])(?:[._-]001)?$", re.IGNORECASE),
    re.compile(r"^(?P<sample>.+?)[._-]?(?:forward|reverse)$", re.IGNORECASE),
]
_FASTQ_ARTIFACTS = (".fastq", ".fq", ".filtered")


def _fastq_stem(name: str) -> str:
    for ext in _FASTQ_EXTS:
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def _read_direction(stem: str) -> int | None:
    """Return 1 for a forward (R1) read, 2 for a reverse (R2) read, else None."""
    for pattern in _READ_PATTERNS:
        match = pattern.match(stem)
        if match is None:
            continue
        read = match.groupdict().get("read")
        if read:
            return int(read)
        return 2 if "reverse" in stem.lower() else 1
    return None


def _first_paired_reads(input_dir: Path) -> tuple[Path | None, Path | None]:
    """Return the first R1 FASTQ and the first R2 FASTQ found in ``input_dir``."""
    r1: Path | None = None
    r2: Path | None = None
    for path in sorted(input_dir.iterdir()):
        if not (path.is_file() and path.name.endswith(_FASTQ_EXTS)):
            continue
        direction = _read_direction(_fastq_stem(path.name))
        if direction == 1 and r1 is None:
            r1 = path
        elif direction == 2 and r2 is None:
            r2 = path
        if r1 is not None and r2 is not None:
            break
    return r1, r2


def _expected_sample_ids(input_dir: Path, *, paired: bool) -> set[str]:
    """Derive expected ASV-table sample ids, mirroring dada2_denoise.R exactly.

    Paired mode strips the R1/R2/read/forward/reverse suffix (R's ``sample_name``).
    Single mode keeps the full stem with no suffix stripping (R's
    ``file_path_sans_ext(file_path_sans_ext(basename(fastqs)))``).
    """
    samples: set[str] = set()
    for path in sorted(input_dir.iterdir()):
        if not (path.is_file() and path.name.endswith(_FASTQ_EXTS)):
            continue
        stem = _fastq_stem(path.name)
        if not paired:
            samples.add(stem)
            continue
        match = next((m for m in (p.match(stem) for p in _READ_PATTERNS) if m), None)
        samples.add(match.group("sample") if match else stem)
    return samples


def _validate_dada2_asv_samples(output_table: Path, input_dir: Path, *, paired: bool) -> None:
    lines = output_table.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise MicrobiomeSuiteError(f"ASV table is empty: {output_table}")
    # write.table(col.names=NA): header is "<empty>\tsample1\tsample2..."
    columns = lines[0].split("\t")[1:]
    if not columns:
        raise MicrobiomeSuiteError(f"ASV table has no sample columns: {output_table}")
    if len(set(columns)) != len(columns):
        raise MicrobiomeSuiteError(f"ASV table has duplicate sample columns: {output_table}")
    for col in columns:
        if not col or any(token in col for token in _FASTQ_ARTIFACTS):
            raise MicrobiomeSuiteError(
                f"ASV table sample column looks like a raw/filtered FASTQ name, not a "
                f"sample id: {col!r} in {output_table}"
            )
    expected = _expected_sample_ids(input_dir, paired=paired)
    actual = set(columns)
    if actual != expected:
        raise MicrobiomeSuiteError(
            f"ASV table sample columns do not match the input samples. "
            f"Unexpected: {sorted(actual - expected)}; missing: {sorted(expected - actual)}."
        )


@dataclass(frozen=True)
class Dada2Tuning:
    """Trim/quality/merge/chimera options shared by the DADA2 backends."""

    trim_left: int = 0
    trunc_len: int = 0
    trim_left_f: int = 0
    trunc_len_f: int = 0
    trim_left_r: int = 0
    trunc_len_r: int = 0
    max_ee: float | None = None
    max_ee_f: float | None = None
    max_ee_r: float | None = None
    trunc_q: int | None = None
    pooling_method: str | None = None
    chimera_method: str | None = None
    min_fold_parent_over_abundance: float | None = None
    allow_one_off: bool | None = None
    n_reads_learn: int | None = None
    homopolymer_gap_penalty: int | None = None
    band_size: int | None = None
    min_overlap: int | None = None
    max_merge_mismatch: int | None = None
    trim_overhang: bool | None = None


def denoise(
    *,
    backend: str,
    demux: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    output_base_transition_stats: Path | None = None,
    output_base_transition_plot: Path | None = None,
    output_plot_dir: Path | None = None,
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
    homopolymer_gap_penalty: int | None = None,
    band_size: int | None = None,
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
    runtime: str = "local",
    dada2_image: str | None = None,
    validate: bool = True,
    amplicon_length: int | None = None,
    strict_qc: bool = False,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "denoise")
    if runtime != "local" and backend != "dada2-r":
        raise MicrobiomeSuiteError("--runtime docker is only supported for --backend dada2-r.")
    if amplicon_length is not None and backend != "dada2-r":
        raise MicrobiomeSuiteError("--amplicon-length only applies to --backend dada2-r.")
    if backend != "dada2-r" and any(
        value is not None for value in (homopolymer_gap_penalty, band_size)
    ):
        raise MicrobiomeSuiteError(
            "--homopolymer-gap-penalty and --band-size only apply to --backend dada2-r."
        )
    resolved_threads = resolve_threads(threads)
    tuning = Dada2Tuning(
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
        homopolymer_gap_penalty=homopolymer_gap_penalty,
        band_size=band_size,
        min_overlap=min_overlap,
        max_merge_mismatch=max_merge_mismatch,
        trim_overhang=trim_overhang,
    )
    if backend == "qiime2-dada2":
        denoise_qiime2_dada2(
            demux=demux,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            output_stats=output_stats,
            output_base_transition_stats=output_base_transition_stats,
            output_base_transition_plot=output_base_transition_plot,
            mode=_resolve_dada2_mode(mode=mode, paired=paired),
            tuning=tuning,
            hashed_feature_ids=hashed_feature_ids,
            retain_all_samples=retain_all_samples,
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
            validate=validate,
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
            validate=validate,
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
            output_plot_dir=output_plot_dir,
            paired=dada2_r_mode == "paired",
            tuning=tuning,
            max_n=max_n,
            rm_phix=rm_phix,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
            runtime=runtime,
            image=dada2_image,
            validate=validate,
            amplicon_length=amplicon_length,
            strict_qc=strict_qc,
        )
        return


def denoise_qiime2_dada2(
    *,
    demux: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    output_base_transition_stats: Path | None,
    output_base_transition_plot: Path | None,
    mode: str,
    tuning: Dada2Tuning,
    hashed_feature_ids: bool | None,
    retain_all_samples: bool | None,
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
    validate: bool = True,
) -> None:
    trim_left = tuning.trim_left
    trunc_len = tuning.trunc_len
    trim_left_f = tuning.trim_left_f
    trunc_len_f = tuning.trunc_len_f
    trim_left_r = tuning.trim_left_r
    trunc_len_r = tuning.trunc_len_r
    max_ee = tuning.max_ee
    max_ee_f = tuning.max_ee_f
    max_ee_r = tuning.max_ee_r
    trunc_q = tuning.trunc_q
    pooling_method = tuning.pooling_method
    chimera_method = tuning.chimera_method
    min_fold_parent_over_abundance = tuning.min_fold_parent_over_abundance
    allow_one_off = tuning.allow_one_off
    n_reads_learn = tuning.n_reads_learn
    min_overlap = tuning.min_overlap
    max_merge_mismatch = tuning.max_merge_mismatch
    trim_overhang = tuning.trim_overhang
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
        inputs={"demux": str(demux)},
        outputs={
            "table": str(output_table),
            "representative_sequences": str(output_rep_seqs),
            "denoising_stats": str(output_stats),
            **(
                {"base_transition_stats": str(output_base_transition_stats)}
                if output_base_transition_stats is not None
                else {}
            ),
        },
        params=params,
        validate=validate,
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
            inputs={"base_transition_stats": str(output_base_transition_stats)},
            outputs={"base_transition_plot": str(output_base_transition_plot)},
            params={**params, "base_transition_plot": str(output_base_transition_plot)},
            validate=validate,
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
    validate: bool = True,
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
        inputs={"demux": str(demux)},
        outputs={
            "table": str(output_table),
            "representative_sequences": str(output_rep_seqs),
            "stats": str(output_stats),
        },
        validate=validate,
    )


def denoise_dada2_r(
    *,
    input_dir: Path,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    output_plot_dir: Path | None,
    paired: bool,
    tuning: Dada2Tuning,
    max_n: int | None,
    rm_phix: bool | None,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
    runtime: str = "local",
    image: str | None = None,
    validate: bool = True,
    amplicon_length: int | None = None,
    strict_qc: bool = False,
) -> None:
    min_overlap = tuning.min_overlap
    max_merge_mismatch = tuning.max_merge_mismatch
    trim_overhang = tuning.trim_overhang
    _validate_dada2_common(
        pooling_method=tuning.pooling_method, chimera_method=tuning.chimera_method
    )
    if not paired and any(
        value is not None for value in (min_overlap, max_merge_mismatch, trim_overhang)
    ):
        raise MicrobiomeSuiteError(
            "--min-overlap, --max-merge-mismatch, and --trim-overhang only apply "
            "to paired DADA2 R mode."
        )

    from importlib.resources import as_file

    from microsuite.runtime.container import (
        PathMapper,
        build_container_command,
        host_user_spec,
        require_engine,
        resolve_dada2_image,
    )

    if runtime == "docker":
        require_engine("docker")
    else:
        rscript = shutil.which("Rscript")
        if rscript is None:
            raise MicrobiomeSuiteError(
                "R/DADA2 denoising requires the external 'Rscript' command. "
                "Install R with the dada2 package, or run it in a container with "
                "--runtime docker (uses the r-dada2 image; no local R needed). "
                "See docs/dada2.md."
            )

    if not input_dir.exists() or not input_dir.is_dir():
        raise MicrobiomeSuiteError(f"Input directory does not exist: {input_dir}")
    _prepare_outputs(output_table, output_rep_seqs, output_stats, force=force)
    if output_plot_dir is not None:
        output_plot_dir.mkdir(parents=True, exist_ok=True)

    overlap_report = None
    if amplicon_length is not None and paired:
        # Inspect one R1 and one R2 independently (mate lengths can differ, e.g.
        # after asymmetric primer trimming), and account for per-mate trim/trunc.
        r1_path, r2_path = _first_paired_reads(input_dir)
        read_len_f = dada2_qc.first_read_length(r1_path) if r1_path is not None else 0
        read_len_r = dada2_qc.first_read_length(r2_path) if r2_path is not None else 0
        overlap_report = dada2_qc.check_overlap(
            trim_left_f=tuning.trim_left_f,
            trim_left_r=tuning.trim_left_r,
            trunc_len_f=tuning.trunc_len_f,
            trunc_len_r=tuning.trunc_len_r,
            read_len_f=read_len_f,
            read_len_r=read_len_r,
            amplicon_length=amplicon_length,
            min_overlap=tuning.min_overlap or 12,
        )
        if overlap_report.warning is not None:
            if strict_qc:
                raise MicrobiomeSuiteError(overlap_report.warning)
            warnings.warn(overlap_report.warning, stacklevel=2)

    script_res = files("microsuite.methods.r").joinpath(DADA2_R_SCRIPT)

    manifest_facts: dict = {
        "microsuite_version": _MICROSUITE_VERSION,
        "backend": "dada2-r",
        "runtime": runtime,
        "image": resolve_dada2_image(image) if runtime == "docker" else None,
        "mode": "paired" if paired else "single",
        "paired": paired,
        "threads": threads,
        "input_dir": str(input_dir),
        "output_table": str(output_table),
        "output_rep_seqs": str(output_rep_seqs),
        "output_stats": str(output_stats),
        "output_plot_dir": str(output_plot_dir) if output_plot_dir is not None else None,
    }
    if overlap_report is not None:
        manifest_facts["overlap_check"] = {
            "read_len_f": overlap_report.read_len_f,
            "read_len_r": overlap_report.read_len_r,
            "retained_f": overlap_report.retained_f,
            "retained_r": overlap_report.retained_r,
            "amplicon_length": overlap_report.amplicon_length,
            "min_overlap": overlap_report.min_overlap,
            "predicted_overlap": overlap_report.predicted_overlap,
            "sufficient": overlap_report.sufficient,
        }

    log_params = _dada2_log_params(
        mode="paired" if paired else "single",
        max_ee=tuning.max_ee,
        max_ee_f=tuning.max_ee_f,
        max_ee_r=tuning.max_ee_r,
        trunc_q=tuning.trunc_q,
        max_n=max_n,
        rm_phix=rm_phix,
        pooling_method=tuning.pooling_method,
        chimera_method=tuning.chimera_method,
        min_fold_parent_over_abundance=tuning.min_fold_parent_over_abundance,
        allow_one_off=tuning.allow_one_off,
        n_reads_learn=tuning.n_reads_learn,
        homopolymer_gap_penalty=tuning.homopolymer_gap_penalty,
        band_size=tuning.band_size,
        min_overlap=min_overlap,
        max_merge_mismatch=max_merge_mismatch,
        trim_overhang=trim_overhang,
    )
    run_kwargs = dict(
        backend="dada2-r",
        inputs={"input_dir": str(input_dir)},
        outputs={
            "table": str(output_table),
            "representative_sequences": str(output_rep_seqs),
            "denoising_stats": str(output_stats),
            **({"plot_dir": str(output_plot_dir)} if output_plot_dir is not None else {}),
        },
        params=log_params,
    )

    stage_dir = output_stats.parent
    with stage_execution(
        stage_dir,
        stage="denoise",
        task="denoise",
        backend="dada2-r",
        params=log_params,
        inputs=[Artifact("reads", input_dir.resolve(), format="directory", required=False)],
    ) as _stage:
        if runtime == "docker":
            # bind mounts cannot expose symlink targets inside input_dir
            for entry in input_dir.iterdir():
                if entry.is_symlink() and entry.name.endswith(
                    (".fastq", ".fq", ".fastq.gz", ".fq.gz")
                ):
                    raise MicrobiomeSuiteError(
                        f"Input FASTQ is a symlink and cannot be mounted into the container: "
                        f"{entry}. Copy/hardlink real files, or use --runtime local."
                    )
            with as_file(script_res) as script_path:
                mapper = PathMapper()
                mapper.add_dir(input_dir, "ro", "/work/input")
                mapper.add_dir(script_path.parent, "ro", "/work/script")
                out_index = 0
                for out in (output_table, output_rep_seqs, output_stats):
                    key = out.resolve().parent
                    if key not in {m.host for m in mapper.mounts()}:
                        mapper.add_dir(key, "rw", f"/work/out{out_index}")
                        out_index += 1
                if output_plot_dir is not None:
                    if output_plot_dir.resolve() not in {m.host for m in mapper.mounts()}:
                        mapper.add_dir(output_plot_dir, "rw", f"/work/out{out_index}")
                inner = [f"{mapper.container_dir(script_path.parent)}/{script_path.name}"]
                inner += _dada2_r_script_args(
                    input_dir=mapper.container_dir(input_dir),
                    output_table=mapper.to_container(output_table),
                    output_rep_seqs=mapper.to_container(output_rep_seqs),
                    output_stats=mapper.to_container(output_stats),
                    output_plot_dir=(
                        mapper.container_dir(output_plot_dir) if output_plot_dir else None
                    ),
                    threads=threads,
                    paired=paired,
                    tuning=tuning,
                    max_n=max_n,
                    rm_phix=rm_phix,
                    params_out=mapper.to_container(output_stats.parent / "dada2_r_params.json"),
                )
                command = build_container_command(
                    inner,
                    resolve_dada2_image(image),
                    mapper.mounts(),
                    engine="docker",
                    user=host_user_spec(),
                )
                _run(
                    command,
                    "R/DADA2 denoising failed.",
                    run_dir=run_dir,
                    timeout=timeout,
                    validate=validate,
                    **run_kwargs,  # ty: ignore[invalid-argument-type]
                )
                if validate:
                    _validate_dada2_asv_samples(output_table, input_dir, paired=paired)
                    summary = dada2_qc.summarize_dada2_stats(output_stats)
                    dada2_qc.write_qc_summary(summary, output_stats.parent)
                    for message in dada2_qc.retention_warnings(summary):
                        if strict_qc:
                            raise MicrobiomeSuiteError(message)
                        warnings.warn(message, stacklevel=2)
                _emit_dada2_manifest(
                    output_stats,
                    {
                        **manifest_facts,
                        "created_at": datetime.now(UTC).isoformat(),
                        "command": " ".join(command),
                    },
                )
        else:
            with as_file(script_res) as script_path:
                command = [rscript, str(script_path)] + _dada2_r_script_args(
                    input_dir=str(input_dir),
                    output_table=str(output_table),
                    output_rep_seqs=str(output_rep_seqs),
                    output_stats=str(output_stats),
                    output_plot_dir=(str(output_plot_dir) if output_plot_dir else None),
                    threads=threads,
                    paired=paired,
                    tuning=tuning,
                    max_n=max_n,
                    rm_phix=rm_phix,
                    params_out=str(output_stats.parent / "dada2_r_params.json"),
                )
                _run(
                    command,
                    "R/DADA2 denoising failed.",
                    run_dir=run_dir,
                    timeout=timeout,
                    validate=validate,
                    **run_kwargs,  # ty: ignore[invalid-argument-type]
                )
                if validate:
                    _validate_dada2_asv_samples(output_table, input_dir, paired=paired)
                    summary = dada2_qc.summarize_dada2_stats(output_stats)
                    dada2_qc.write_qc_summary(summary, output_stats.parent)
                    for message in dada2_qc.retention_warnings(summary):
                        if strict_qc:
                            raise MicrobiomeSuiteError(message)
                        warnings.warn(message, stacklevel=2)
                _emit_dada2_manifest(
                    output_stats,
                    {
                        **manifest_facts,
                        "created_at": datetime.now(UTC).isoformat(),
                        "command": " ".join(command),
                    },
                )
        _declare_dada2_stage_outputs(
            _stage,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            output_stats=output_stats,
            output_plot_dir=output_plot_dir,
            manifest_path=stage_dir / dada2_manifest.MANIFEST_FILENAME,
            validate=validate,
        )


def _declare_dada2_stage_outputs(
    stage: StageRecord,
    *,
    output_table: Path,
    output_rep_seqs: Path,
    output_stats: Path,
    output_plot_dir: Path | None,
    manifest_path: Path,
    validate: bool,
) -> None:
    """Declare the dada2-r outputs + provenance on the active stage record.

    Required-ness of the primary outputs mirrors ``validate``: a caller that
    skipped output validation is not asserting the files exist. The provenance
    manifest is best-effort (``required=False``) because ``_emit_dada2_manifest``
    tolerates a missing DADA2 params file.
    """
    stage.add_output(
        Artifact(
            "feature table",
            output_table.resolve(),
            format="tsv",
            kind="feature_table",
            required=validate,
        )
    )
    stage.add_output(
        Artifact(
            "representative sequences",
            output_rep_seqs.resolve(),
            format="fasta",
            kind="representative_sequences",
            required=validate,
        )
    )
    stage.add_output(
        Artifact(
            "denoising stats",
            output_stats.resolve(),
            format="tsv",
            kind="denoising_stats",
            required=validate,
        )
    )
    if output_plot_dir is not None:
        stage.add_output(
            Artifact(
                "plots",
                output_plot_dir.resolve(),
                format="directory",
                kind="plots",
                required=False,
            )
        )
    stage.add_provenance(ProvenanceFile("dada2_manifest", manifest_path.resolve(), required=False))
    _lift_dada2_software(stage, manifest_path)
    _lift_dada2_metrics(stage, output_stats)


def _lift_dada2_metrics(stage: StageRecord, output_stats: Path) -> None:
    """Best-effort: emit read-retention metrics from the DADA2 stats summary."""
    if not output_stats.exists():
        return
    try:
        overall = dada2_qc.summarize_dada2_stats(output_stats)["overall"]
    except (MicrobiomeSuiteError, OSError, KeyError, ValueError):
        return
    inp = overall.get("input") or 0
    stage.add_metric("input_reads", inp, "reads")
    if overall.get("filtered_frac") is not None:
        stage.add_metric("filtered_fraction", overall["filtered_frac"], "fraction")
    if overall.get("merged_frac") is not None:
        stage.add_metric("merged_fraction", overall["merged_frac"], "fraction")
    nonchim = overall.get("nonchim_frac")
    if nonchim is not None:
        stage.add_metric("nonchimeric_fraction", nonchim, "fraction")
        stage.add_metric("nonchimeric_reads", round(inp * nonchim), "reads")


def _lift_dada2_software(stage: StageRecord, manifest_path: Path) -> None:
    """Best-effort: declare dada2/R/microsuite versions from the manifest tool block.

    Reads the persisted ``dada2_denoise_manifest.json`` (``_emit_dada2_manifest``
    deletes the intermediate R params file), so a missing/unparseable manifest
    simply leaves ``software`` empty and never fails the run.
    """
    if not manifest_path.exists():
        return
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return
    tool = manifest.get("tool", {}) if isinstance(manifest, dict) else {}
    software: dict[str, dict[str, str]] = {"microsuite": {"version": _MICROSUITE_VERSION}}
    for key, name in (("dada2_version", "dada2"), ("r_version", "R")):
        version = tool.get(key)
        if version is not None:
            software[name] = {"version": str(version)}
    stage.set_software(software)


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


def _validate_dada2_common(*, pooling_method: str | None, chimera_method: str | None) -> None:
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


def _dada2_r_script_args(
    *,
    input_dir: str,
    output_table: str,
    output_rep_seqs: str,
    output_stats: str,
    output_plot_dir: str | None,
    threads: int,
    paired: bool,
    tuning: Dada2Tuning,
    max_n: int | None,
    rm_phix: bool | None,
    params_out: str | None = None,
) -> list[str]:
    args = [
        "--input-dir",
        input_dir,
        "--output-table",
        output_table,
        "--output-rep-seqs",
        output_rep_seqs,
        "--output-stats",
        output_stats,
        "--threads",
        str(threads),
    ]
    if output_plot_dir is not None:
        args.extend(["--output-plot-dir", output_plot_dir])
    if paired:
        args.append("--paired")
        args.extend(
            [
                "--trim-left-f",
                str(tuning.trim_left_f),
                "--trunc-len-f",
                str(tuning.trunc_len_f),
                "--trim-left-r",
                str(tuning.trim_left_r),
                "--trunc-len-r",
                str(tuning.trunc_len_r),
            ]
        )
        _append_value(args, "--max-ee-f", tuning.max_ee_f)
        _append_value(args, "--max-ee-r", tuning.max_ee_r)
        _append_value(args, "--min-overlap", tuning.min_overlap)
        _append_value(args, "--max-merge-mismatch", tuning.max_merge_mismatch)
        _append_bool(args, "--trim-overhang", tuning.trim_overhang)
    else:
        args.extend(["--trim-left", str(tuning.trim_left), "--trunc-len", str(tuning.trunc_len)])
        _append_value(args, "--max-ee", tuning.max_ee)
    _append_value(args, "--trunc-q", tuning.trunc_q)
    _append_value(args, "--max-n", max_n)
    _append_bool(args, "--rm-phix", rm_phix)
    _append_value(args, "--pooling-method", tuning.pooling_method)
    _append_value(args, "--chimera-method", tuning.chimera_method)
    _append_value(args, "--min-fold-parent-over-abundance", tuning.min_fold_parent_over_abundance)
    _append_bool(args, "--allow-one-off", tuning.allow_one_off)
    _append_value(args, "--n-reads-learn", tuning.n_reads_learn)
    _append_value(args, "--homopolymer-gap-penalty", tuning.homopolymer_gap_penalty)
    _append_value(args, "--band-size", tuning.band_size)
    if params_out is not None:
        args.extend(["--params-out", params_out])
    return args


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


def _emit_dada2_manifest(output_stats: Path, run_facts: dict) -> None:
    """Merge the R-emitted resolved params with wrapper facts into
    dada2_denoise_manifest.json beside the stats file. Best-effort: a missing or
    unparseable R params file warns and never fails an otherwise-successful run."""
    r_params_path = output_stats.parent / "dada2_r_params.json"
    try:
        r_params = dada2_manifest.read_r_params(r_params_path)
        manifest = dada2_manifest.build_manifest(r_params, run_facts)
        dada2_manifest.write_manifest(manifest, output_stats.parent)
    except MicrobiomeSuiteError as exc:
        warnings.warn(f"Could not write DADA2 provenance manifest: {exc}", stacklevel=2)
    finally:
        r_params_path.unlink(missing_ok=True)


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
    inputs: dict[str, str] | None = None,
    outputs: dict[str, str] | None = None,
    params: dict[str, object] | None = None,
    validate: bool = True,
) -> None:
    run_command(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="denoise",
            backend=backend,
            inputs=inputs,
            outputs=outputs,
            params=params or {},
        ),
    )
    if validate and outputs:
        validate_outputs(outputs)
