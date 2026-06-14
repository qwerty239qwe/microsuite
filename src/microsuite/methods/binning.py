from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("metabat2", "maxbin2", "concoct")


def bin_contigs(
    *,
    backend: str,
    contigs: Path,
    output_dir: Path,
    depth: Path | None = None,
    abundance: Path | None = None,
    coverage: Path | None = None,
    prefix: str = "bin",
    threads: int | str = "1",
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "bin")
    if backend == "metabat2":
        bin_metabat2(
            contigs=contigs,
            depth=depth,
            output_dir=output_dir,
            prefix=prefix,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "maxbin2":
        bin_maxbin2(
            contigs=contigs,
            abundance=abundance,
            output_dir=output_dir,
            prefix=prefix,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "concoct":
        bin_concoct(
            contigs=contigs,
            coverage=coverage,
            output_dir=output_dir,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return


def bin_metabat2(
    *,
    contigs: Path,
    depth: Path | None,
    output_dir: Path,
    prefix: str,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if depth is None:
        raise MicrobiomeSuiteError("--depth is required for --backend metabat2.")
    executable = _require_tool(
        "metabat2",
        "MetaBAT2 binning requires the external 'metabat2' command.",
    )
    output_dir = _prepare_output_dir(output_dir, force=force)
    output_prefix = output_dir / prefix
    command = [
        executable,
        "-i",
        str(ensure_input(contigs)),
        "-a",
        str(ensure_input(depth)),
        "-o",
        str(output_prefix),
        "-t",
        str(resolve_threads(threads)),
    ]
    _run(
        command,
        "MetaBAT2 binning failed.",
        backend="metabat2",
        inputs={"contigs": str(contigs), "depth": str(depth)},
        output_dir=output_dir,
        params={"prefix": prefix},
        run_dir=run_dir,
        timeout=timeout,
    )


def bin_maxbin2(
    *,
    contigs: Path,
    abundance: Path | None,
    output_dir: Path,
    prefix: str,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if abundance is None:
        raise MicrobiomeSuiteError("--abundance is required for --backend maxbin2.")
    executable = _require_tool(
        "run_MaxBin.pl",
        "MaxBin2 binning requires the external 'run_MaxBin.pl' command.",
    )
    output_dir = _prepare_output_dir(output_dir, force=force)
    output_prefix = output_dir / prefix
    command = [
        executable,
        "-contig",
        str(ensure_input(contigs)),
        "-abund",
        str(ensure_input(abundance)),
        "-out",
        str(output_prefix),
        "-thread",
        str(resolve_threads(threads)),
    ]
    _run(
        command,
        "MaxBin2 binning failed.",
        backend="maxbin2",
        inputs={"contigs": str(contigs), "abundance": str(abundance)},
        output_dir=output_dir,
        params={"prefix": prefix},
        run_dir=run_dir,
        timeout=timeout,
    )


def bin_concoct(
    *,
    contigs: Path,
    coverage: Path | None,
    output_dir: Path,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if coverage is None:
        raise MicrobiomeSuiteError("--coverage is required for --backend concoct.")
    executable = _require_tool(
        "concoct",
        "CONCOCT binning requires the external 'concoct' command.",
    )
    output_dir = _prepare_output_dir(output_dir, force=force)
    command = [
        executable,
        "--composition_file",
        str(ensure_input(contigs)),
        "--coverage_file",
        str(ensure_input(coverage)),
        "-b",
        str(output_dir),
        "-t",
        str(resolve_threads(threads)),
    ]
    _run(
        command,
        "CONCOCT binning failed.",
        backend="concoct",
        inputs={"contigs": str(contigs), "coverage": str(coverage)},
        output_dir=output_dir,
        params={},
        run_dir=run_dir,
        timeout=timeout,
    )


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
    params: dict[str, str],
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    run_command(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="bin",
            backend=backend,
            inputs=inputs,
            outputs={"output_dir": str(output_dir)},
            params=params,
        ),
    )
