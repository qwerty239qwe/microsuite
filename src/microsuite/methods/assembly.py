from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("megahit", "metaspades", "idba-ud", "mosh-megahit")


def assemble(
    *,
    backend: str,
    output_dir: Path | None = None,
    read1: Path | None = None,
    read2: Path | None = None,
    reads: Path | None = None,
    output_contigs: Path | None = None,
    presets: str = "meta-sensitive",
    min_contig: int = 500,
    parallel_config: Path | None = None,
    verbose: bool = False,
    threads: int | str = "1",
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "assemble")
    if backend == "megahit":
        assemble_megahit(
            read1=read1,
            read2=read2,
            reads=reads,
            output_dir=output_dir,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "metaspades":
        assemble_metaspades(
            read1=read1,
            read2=read2,
            reads=reads,
            output_dir=output_dir,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "idba-ud":
        assemble_idba_ud(
            reads=reads,
            output_dir=output_dir,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "mosh-megahit":
        assemble_mosh_megahit(
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
        return


def assemble_megahit(
    *,
    read1: Path | None,
    read2: Path | None,
    reads: Path | None,
    output_dir: Path,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    output_dir = _require_output_dir(output_dir)
    _validate_assembly_reads(read1=read1, read2=read2, reads=reads)
    executable = _require_tool(
        "megahit",
        "MEGAHIT assembly requires the external 'megahit' command.",
    )
    output_dir = _prepare_output_dir(output_dir, force=force)
    command = [executable]
    inputs = _append_assembly_reads(
        command,
        read1=read1,
        read2=read2,
        reads=reads,
        single_flag="-r",
    )
    command.extend(["-o", str(output_dir), "-t", str(resolve_threads(threads))])
    _run(
        command,
        "MEGAHIT assembly failed.",
        backend="megahit",
        inputs=inputs,
        output_dir=output_dir,
        run_dir=run_dir,
        timeout=timeout,
    )


def assemble_metaspades(
    *,
    read1: Path | None,
    read2: Path | None,
    reads: Path | None,
    output_dir: Path,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    output_dir = _require_output_dir(output_dir)
    _validate_assembly_reads(read1=read1, read2=read2, reads=reads)
    executable = _require_tool(
        "metaspades.py",
        "metaSPAdes assembly requires the external 'metaspades.py' command.",
    )
    output_dir = _prepare_output_dir(output_dir, force=force)
    command = [executable]
    inputs = _append_assembly_reads(
        command,
        read1=read1,
        read2=read2,
        reads=reads,
        single_flag="-s",
    )
    command.extend(["-o", str(output_dir), "-t", str(resolve_threads(threads))])
    _run(
        command,
        "metaSPAdes assembly failed.",
        backend="metaspades",
        inputs=inputs,
        output_dir=output_dir,
        run_dir=run_dir,
        timeout=timeout,
    )


def assemble_idba_ud(
    *,
    reads: Path | None,
    output_dir: Path,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    output_dir = _require_output_dir(output_dir)
    if reads is None:
        raise MicrobiomeSuiteError("--reads is required for --backend idba-ud.")
    executable = _require_tool(
        "idba_ud",
        "IDBA-UD assembly requires the external 'idba_ud' command.",
    )
    output_dir = _prepare_output_dir(output_dir, force=force)
    command = [
        executable,
        "-r",
        str(ensure_input(reads)),
        "-o",
        str(output_dir),
        "--num_threads",
        str(resolve_threads(threads)),
    ]
    _run(
        command,
        "IDBA-UD assembly failed.",
        backend="idba-ud",
        inputs={"reads": str(reads)},
        output_dir=output_dir,
        run_dir=run_dir,
        timeout=timeout,
    )


def assemble_mosh_megahit(
    *,
    reads: Path | None,
    output_dir: Path | None,
    output_contigs: Path | None,
    presets: str,
    min_contig: int,
    parallel_config: Path | None,
    verbose: bool,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if reads is None:
        raise MicrobiomeSuiteError("--reads is required for --backend mosh-megahit.")
    if min_contig < 1:
        raise MicrobiomeSuiteError("--min-contig must be greater than zero.")
    executable = _require_tool(
        "mosh",
        "MOSHPIT assembly requires the external 'mosh' command.",
    )
    output_contigs = _resolve_output_artifact(
        output=output_contigs,
        output_dir=output_dir,
        default_name="contigs.qza",
        force=force,
    )
    command = [
        executable,
        "assembly",
        "assemble-megahit",
        "--i-reads",
        str(ensure_input(reads)),
        "--p-presets",
        presets,
        "--p-num-cpu-threads",
        str(resolve_threads(threads)),
        "--p-min-contig",
        str(min_contig),
        "--o-contigs",
        str(output_contigs),
    ]
    if parallel_config is not None:
        command.extend(["--parallel-config", str(ensure_input(parallel_config))])
    if verbose:
        command.append("--verbose")
    _run(
        command,
        "MOSHPIT MEGAHIT assembly failed.",
        backend="mosh-megahit",
        inputs={"reads": str(reads)},
        outputs={"contigs": str(output_contigs)},
        params={
            "presets": presets,
            "min_contig": str(min_contig),
            "threads": str(resolve_threads(threads)),
        },
        run_dir=run_dir,
        timeout=timeout,
    )


def _validate_assembly_reads(
    *,
    read1: Path | None,
    read2: Path | None,
    reads: Path | None,
) -> None:
    if reads is not None and (read1 is not None or read2 is not None):
        raise MicrobiomeSuiteError("Use either --reads or --read1/--read2, not both.")
    if reads is None and read1 is None:
        raise MicrobiomeSuiteError("--read1 or --reads is required for assembly.")
    if read1 is None and read2 is not None:
        raise MicrobiomeSuiteError("--read2 requires --read1.")


def _append_assembly_reads(
    command: list[str],
    *,
    read1: Path | None,
    read2: Path | None,
    reads: Path | None,
    single_flag: str,
) -> dict[str, str]:
    if reads is not None:
        command.extend([single_flag, str(ensure_input(reads))])
        return {"reads": str(reads)}
    if read1 is None:
        raise MicrobiomeSuiteError("--read1 or --reads is required for assembly.")
    if read2 is None:
        command.extend([single_flag, str(ensure_input(read1))])
        return {"read1": str(read1)}
    command.extend(["-1", str(ensure_input(read1)), "-2", str(ensure_input(read2))])
    return {"read1": str(read1), "read2": str(read2)}


def _require_tool(name: str, message: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise MicrobiomeSuiteError(message)
    return executable


def _require_output_dir(path: Path | None) -> Path:
    if path is None:
        raise MicrobiomeSuiteError("--output-dir is required for this assembly backend.")
    return path


def _prepare_output_dir(path: Path, *, force: bool) -> Path:
    if path.exists():
        if not path.is_dir():
            raise MicrobiomeSuiteError(f"Output path exists and is not a directory: {path}")
        if any(path.iterdir()) and not force:
            raise MicrobiomeSuiteError(
                f"Output directory exists, pass --force to overwrite: {path}"
            )
    path.mkdir(parents=True, exist_ok=True)
    return path


def _resolve_output_artifact(
    *,
    output: Path | None,
    output_dir: Path | None,
    default_name: str,
    force: bool,
) -> Path:
    if output is None:
        if output_dir is None:
            raise MicrobiomeSuiteError(
                f"--output-contigs or --output-dir is required to write {default_name}."
            )
        output_dir = _prepare_output_dir(output_dir, force=force)
        output = output_dir / default_name
    elif output.exists() and not force:
        raise MicrobiomeSuiteError(f"Output file exists, pass --force to overwrite: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def _run(
    command: list[str],
    failure_message: str,
    *,
    backend: str,
    inputs: dict[str, str],
    output_dir: Path | None = None,
    outputs: dict[str, str] | None = None,
    params: dict[str, str] | None = None,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    logged_outputs = outputs or {"output_dir": str(output_dir)}
    run_command(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="assemble",
            backend=backend,
            inputs=inputs,
            outputs=logged_outputs,
            params=params,
        ),
    )
