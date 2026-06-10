from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("megahit", "metaspades", "idba-ud")


def assemble(
    *,
    backend: str,
    output_dir: Path,
    read1: Path | None = None,
    read2: Path | None = None,
    reads: Path | None = None,
    threads: int | str = "1",
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
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
    raise MicrobiomeSuiteError(
        f"Unsupported assemble backend '{backend}'. Choose one of: {', '.join(SUPPORTED_BACKENDS)}"
    )


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


def _run(
    command: list[str],
    failure_message: str,
    *,
    backend: str,
    inputs: dict[str, str],
    output_dir: Path,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    run_command(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="assemble",
            backend=backend,
            inputs=inputs,
            outputs={"output_dir": str(output_dir)},
        ),
    )
