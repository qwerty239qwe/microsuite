from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output

SUPPORTED_BACKENDS = ("fastp", "cutadapt", "qiime2-cutadapt")
PLANNED_BACKENDS = ("cutadapt", "qiime2-cutadapt")


def trim(
    *,
    backend: str,
    read1: Path,
    output1: Path,
    read2: Path | None = None,
    output2: Path | None = None,
    html: Path | None = None,
    json_report: Path | None = None,
    threads: int = 1,
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend == "fastp":
        trim_fastp(
            read1=read1,
            output1=output1,
            read2=read2,
            output2=output2,
            html=html,
            json_report=json_report,
            threads=threads,
            force=force,
        )
        return
    if backend in PLANNED_BACKENDS:
        raise MicrobiomeSuiteError(
            f"Trim backend '{backend}' is registered but not implemented yet. "
            "Use --backend fastp for now."
        )
    backends = ", ".join(SUPPORTED_BACKENDS)
    raise MicrobiomeSuiteError(f"Unsupported trim backend '{backend}'. Choose one of: {backends}")


def trim_fastp(
    *,
    read1: Path,
    output1: Path,
    read2: Path | None,
    output2: Path | None,
    html: Path | None,
    json_report: Path | None,
    threads: int,
    force: bool,
) -> None:
    if read2 is not None and output2 is None:
        raise MicrobiomeSuiteError("--output2 is required when --read2 is supplied.")
    if output2 is not None and read2 is None:
        raise MicrobiomeSuiteError("--read2 is required when --output2 is supplied.")
    fastp = shutil.which("fastp")
    if fastp is None:
        raise MicrobiomeSuiteError(
            "fastp trimming requires the external 'fastp' command. "
            "Install fastp and rerun this command."
        )

    ensure_input(read1)
    if read2 is not None:
        ensure_input(read2)
    _prepare_outputs(output1, output2, html, json_report, force=force)

    command = [fastp, "--in1", str(read1), "--out1", str(output1)]
    if read2 is not None and output2 is not None:
        command.extend(["--in2", str(read2), "--out2", str(output2)])
    if html is not None:
        command.extend(["--html", str(html)])
    if json_report is not None:
        command.extend(["--json", str(json_report)])
    command.extend(["--thread", str(threads)])

    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or "fastp trimming failed."
        raise MicrobiomeSuiteError(message)


def _prepare_outputs(*outputs: Path | None, force: bool) -> None:
    for output in outputs:
        if output is not None:
            prepare_output(output, force=force)
