from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("metabat2", "maxbin2", "concoct", "mosh-metabat2")


def bin_contigs(
    *,
    backend: str,
    contigs: Path,
    output_dir: Path | None = None,
    depth: Path | None = None,
    abundance: Path | None = None,
    coverage: Path | None = None,
    alignment_maps: Path | None = None,
    output_mags: Path | None = None,
    output_contig_map: Path | None = None,
    output_unbinned_contigs: Path | None = None,
    prefix: str = "bin",
    seed: int = 100,
    parallel_config: Path | None = None,
    verbose: bool = False,
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
            output_dir=_require_output_dir(output_dir),
            prefix=prefix,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "mosh-metabat2":
        bin_mosh_metabat2(
            contigs=contigs,
            alignment_maps=alignment_maps,
            output_dir=output_dir,
            output_mags=output_mags,
            output_contig_map=output_contig_map,
            output_unbinned_contigs=output_unbinned_contigs,
            seed=seed,
            parallel_config=parallel_config,
            verbose=verbose,
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
            output_dir=_require_output_dir(output_dir),
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
            output_dir=_require_output_dir(output_dir),
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
    output_dir = _require_output_dir(output_dir)
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
    output_dir = _require_output_dir(output_dir)
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
    output_dir = _require_output_dir(output_dir)
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


def bin_mosh_metabat2(
    *,
    contigs: Path,
    alignment_maps: Path | None,
    output_dir: Path | None,
    output_mags: Path | None,
    output_contig_map: Path | None,
    output_unbinned_contigs: Path | None,
    seed: int,
    parallel_config: Path | None,
    verbose: bool,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if alignment_maps is None:
        raise MicrobiomeSuiteError("--alignment-maps is required for --backend mosh-metabat2.")
    if seed < 0:
        raise MicrobiomeSuiteError("--seed must be greater than or equal to zero.")
    executable = _require_tool(
        "mosh",
        "MOSHPIT binning requires the external 'mosh' command.",
    )
    outputs = _resolve_mosh_outputs(
        output_dir=output_dir,
        output_mags=output_mags,
        output_contig_map=output_contig_map,
        output_unbinned_contigs=output_unbinned_contigs,
        force=force,
    )
    command = [
        executable,
        "annotate",
        "bin-contigs-metabat",
        "--i-contigs",
        str(ensure_input(contigs)),
        "--i-alignment-maps",
        str(ensure_input(alignment_maps)),
        "--p-num-threads",
        str(resolve_threads(threads)),
        "--p-seed",
        str(seed),
        "--o-mags",
        str(outputs["mags"]),
        "--o-contig-map",
        str(outputs["contig_map"]),
        "--o-unbinned-contigs",
        str(outputs["unbinned_contigs"]),
    ]
    if parallel_config is not None:
        command.extend(["--parallel-config", str(ensure_input(parallel_config))])
    if verbose:
        command.append("--verbose")
    _run(
        command,
        "MOSHPIT MetaBAT2 binning failed.",
        backend="mosh-metabat2",
        inputs={"contigs": str(contigs), "alignment_maps": str(alignment_maps)},
        outputs={key: str(value) for key, value in outputs.items()},
        params={"seed": str(seed), "threads": str(resolve_threads(threads))},
        run_dir=run_dir,
        timeout=timeout,
    )


def _require_tool(name: str, message: str) -> str:
    executable = shutil.which(name)
    if executable is None:
        raise MicrobiomeSuiteError(message)
    return executable


def _require_output_dir(path: Path | None) -> Path:
    if path is None:
        raise MicrobiomeSuiteError("--output-dir is required for this binning backend.")
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


def _resolve_mosh_outputs(
    *,
    output_dir: Path | None,
    output_mags: Path | None,
    output_contig_map: Path | None,
    output_unbinned_contigs: Path | None,
    force: bool,
) -> dict[str, Path]:
    if output_dir is None and any(
        output is None for output in (output_mags, output_contig_map, output_unbinned_contigs)
    ):
        raise MicrobiomeSuiteError(
            "--output-dir or all MOSHPIT output artifacts are required for mosh-metabat2."
        )
    resolved_output_dir = (
        _prepare_output_dir(output_dir, force=force) if output_dir is not None else None
    )
    outputs = {
        "mags": output_mags or _join_output_dir(resolved_output_dir, "mags.qza"),
        "contig_map": output_contig_map or _join_output_dir(resolved_output_dir, "contig-map.qza"),
        "unbinned_contigs": output_unbinned_contigs
        or _join_output_dir(resolved_output_dir, "unbinned-contigs.qza"),
    }
    for output in outputs.values():
        if output.exists() and not force:
            raise MicrobiomeSuiteError(f"Output file exists, pass --force to overwrite: {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
    return outputs


def _join_output_dir(output_dir: Path | None, filename: str) -> Path:
    if output_dir is None:
        raise MicrobiomeSuiteError(
            "--output-dir or all MOSHPIT output artifacts are required for mosh-metabat2."
        )
    return output_dir / filename


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
            task="bin",
            backend=backend,
            inputs=inputs,
            outputs=logged_outputs,
            params=params,
        ),
    )
