from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output

SUPPORTED_BACKENDS = ("fastqc", "multiqc", "qiime2-demux")


def qc(
    *,
    backend: str,
    inputs: list[Path] | None = None,
    input_dir: Path | None = None,
    demux: Path | None = None,
    output_dir: Path | None = None,
    output: Path | None = None,
    threads: int = 1,
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend == "fastqc":
        qc_fastqc(inputs=inputs or [], output_dir=output_dir, threads=threads, force=force)
        return
    if backend == "multiqc":
        qc_multiqc(input_dir=input_dir, output_dir=output_dir, force=force)
        return
    if backend == "qiime2-demux":
        qc_qiime2_demux(demux=demux, output=output, force=force)
        return
    backends = ", ".join(SUPPORTED_BACKENDS)
    raise MicrobiomeSuiteError(f"Unsupported QC backend '{backend}'. Choose one of: {backends}")


def qc_fastqc(*, inputs: list[Path], output_dir: Path | None, threads: int, force: bool) -> None:
    if not inputs:
        raise MicrobiomeSuiteError("--input is required for --backend fastqc.")
    if output_dir is None:
        raise MicrobiomeSuiteError("--output-dir is required for --backend fastqc.")
    fastqc = _require_tool("fastqc", "FastQC")
    for path in inputs:
        ensure_input(path)
    _prepare_dir(output_dir, force=force)

    command = [fastqc, "--outdir", str(output_dir), "--threads", str(threads)]
    command.extend(str(path) for path in inputs)
    _run(command, "FastQC failed.")


def qc_multiqc(*, input_dir: Path | None, output_dir: Path | None, force: bool) -> None:
    if input_dir is None:
        raise MicrobiomeSuiteError("--input-dir is required for --backend multiqc.")
    if output_dir is None:
        raise MicrobiomeSuiteError("--output-dir is required for --backend multiqc.")
    multiqc = _require_tool("multiqc", "MultiQC")
    if not input_dir.exists() or not input_dir.is_dir():
        raise MicrobiomeSuiteError(f"Input directory does not exist: {input_dir}")
    _prepare_dir(output_dir, force=force)

    command = [multiqc, str(input_dir), "--outdir", str(output_dir)]
    if force:
        command.append("--force")
    _run(command, "MultiQC failed.")


def qc_qiime2_demux(*, demux: Path | None, output: Path | None, force: bool) -> None:
    if demux is None:
        raise MicrobiomeSuiteError("--demux is required for --backend qiime2-demux.")
    if output is None:
        raise MicrobiomeSuiteError("--output is required for --backend qiime2-demux.")
    qiime = _require_tool("qiime", "QIIME 2 demux summary")
    ensure_input(demux)
    if output.exists() and force:
        output.unlink()
    prepare_output(output, force=force)

    command = [
        qiime,
        "demux",
        "summarize",
        "--i-data",
        str(demux),
        "--o-visualization",
        str(output),
    ]
    _run(command, "QIIME 2 demux summary failed.")


def _require_tool(command: str, task: str) -> str:
    executable = shutil.which(command)
    if executable is None:
        raise MicrobiomeSuiteError(
            f"{task} requires the external '{command}' command. "
            "Install or activate the tool and rerun this command."
        )
    return executable


def _prepare_dir(path: Path, *, force: bool) -> None:
    if path.exists() and not path.is_dir():
        raise MicrobiomeSuiteError(f"Output path exists and is not a directory: {path}")
    if path.exists() and any(path.iterdir()) and not force:
        raise MicrobiomeSuiteError(f"Output directory is not empty, pass --force to use it: {path}")
    path.mkdir(parents=True, exist_ok=True)


def _run(command: list[str], failure_message: str) -> None:
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or failure_message
        raise MicrobiomeSuiteError(message)
