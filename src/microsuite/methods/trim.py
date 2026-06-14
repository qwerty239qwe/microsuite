from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.methods._dispatch import require_backend
from microsuite.methods._qiime import require_qiime, run_qiime
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("fastp", "cutadapt", "trimmomatic", "trim-galore", "qiime2-cutadapt")


def trim(
    *,
    backend: str,
    read1: Path,
    output1: Path,
    read2: Path | None = None,
    output2: Path | None = None,
    unpaired1: Path | None = None,
    unpaired2: Path | None = None,
    html: Path | None = None,
    json_report: Path | None = None,
    adapter: str | None = None,
    front: str | None = None,
    anywhere: str | None = None,
    adapter2: str | None = None,
    front2: str | None = None,
    anywhere2: str | None = None,
    quality_cutoff: str | None = None,
    minimum_length: str | None = None,
    maximum_length: str | None = None,
    max_n: str | None = None,
    discard_untrimmed: bool = False,
    trimmomatic_steps: list[str] | None = None,
    basename: str | None = None,
    trim_galore_version: str = "auto",
    threads: int | str = 1,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "trim")
    resolved_threads = resolve_threads(threads)
    if backend == "fastp":
        if (
            any(
                value is not None
                for value in (
                    front,
                    anywhere,
                    front2,
                    anywhere2,
                )
            )
            or discard_untrimmed
        ):
            raise MicrobiomeSuiteError("Cutadapt-specific trim options require --backend cutadapt.")
        trim_fastp(
            read1=read1,
            output1=output1,
            read2=read2,
            output2=output2,
            html=html,
            json_report=json_report,
            adapter=adapter,
            adapter2=adapter2,
            quality_cutoff=quality_cutoff,
            minimum_length=minimum_length,
            maximum_length=maximum_length,
            max_n=max_n,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "cutadapt":
        trim_cutadapt(
            read1=read1,
            output1=output1,
            read2=read2,
            output2=output2,
            unpaired1=unpaired1,
            unpaired2=unpaired2,
            html=html,
            json_report=json_report,
            adapter=adapter,
            front=front,
            anywhere=anywhere,
            adapter2=adapter2,
            front2=front2,
            anywhere2=anywhere2,
            quality_cutoff=quality_cutoff,
            minimum_length=minimum_length,
            maximum_length=maximum_length,
            max_n=max_n,
            discard_untrimmed=discard_untrimmed,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "trimmomatic":
        _reject_options(
            "trimmomatic",
            {
                "--html": html,
                "--json-report": json_report,
                "--adapter": adapter,
                "--front": front,
                "--anywhere": anywhere,
                "--adapter2": adapter2,
                "--front2": front2,
                "--anywhere2": anywhere2,
                "--quality-cutoff": quality_cutoff,
                "--minimum-length": minimum_length,
                "--maximum-length": maximum_length,
                "--max-n": max_n,
                "--discard-untrimmed": discard_untrimmed,
                "--basename": basename,
            },
        )
        trim_trimmomatic(
            read1=read1,
            output1=output1,
            read2=read2,
            output2=output2,
            unpaired1=unpaired1,
            unpaired2=unpaired2,
            steps=trimmomatic_steps or [],
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "trim-galore":
        _reject_options(
            "trim-galore",
            {
                "--html": html,
                "--json-report": json_report,
                "--front": front,
                "--anywhere": anywhere,
                "--front2": front2,
                "--anywhere2": anywhere2,
                "--maximum-length": maximum_length,
                "--max-n": max_n,
                "--discard-untrimmed": discard_untrimmed,
                "--trimmomatic-step": trimmomatic_steps,
                "--unpaired1": unpaired1,
                "--unpaired2": unpaired2,
            },
        )
        trim_galore(
            read1=read1,
            output1=output1,
            read2=read2,
            output2=output2,
            adapter=adapter,
            adapter2=adapter2,
            quality_cutoff=quality_cutoff,
            minimum_length=minimum_length,
            basename=basename,
            trim_galore_version=trim_galore_version,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "qiime2-cutadapt":
        _reject_options(
            "qiime2-cutadapt",
            {
                "--read2": read2,
                "--output2": output2,
                "--unpaired1": unpaired1,
                "--unpaired2": unpaired2,
                "--html": html,
                "--json-report": json_report,
                "--basename": basename,
                "--trimmomatic-step": trimmomatic_steps,
            },
        )
        trim_qiime2_cutadapt(
            demux=read1,
            output=output1,
            adapter=adapter,
            front=front,
            anywhere=anywhere,
            adapter2=adapter2,
            front2=front2,
            anywhere2=anywhere2,
            quality_cutoff=quality_cutoff,
            minimum_length=minimum_length,
            maximum_length=maximum_length,
            max_n=max_n,
            discard_untrimmed=discard_untrimmed,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return


def trim_fastp(
    *,
    read1: Path,
    output1: Path,
    read2: Path | None,
    output2: Path | None,
    html: Path | None,
    json_report: Path | None,
    adapter: str | None,
    adapter2: str | None,
    quality_cutoff: str | None,
    minimum_length: str | None,
    maximum_length: str | None,
    max_n: str | None,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if read2 is not None and output2 is None:
        raise MicrobiomeSuiteError("--output2 is required when --read2 is supplied.")
    if output2 is not None and read2 is None:
        raise MicrobiomeSuiteError("--read2 is required when --output2 is supplied.")
    if adapter2 is not None and read2 is None:
        raise MicrobiomeSuiteError("--adapter2 requires --read2 and --output2 for --backend fastp.")
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
    if adapter is not None:
        command.extend(["--adapter_sequence", adapter])
    if adapter2 is not None:
        command.extend(["--adapter_sequence_r2", adapter2])
    _add_optional(command, "--qualified_quality_phred", quality_cutoff)
    _add_optional(command, "--length_required", minimum_length)
    _add_optional(command, "--length_limit", maximum_length)
    _add_optional(command, "--n_base_limit", max_n)
    if html is not None:
        command.extend(["--html", str(html)])
    if json_report is not None:
        command.extend(["--json", str(json_report)])
    command.extend(["--thread", str(threads)])

    _run(command, "fastp trimming failed.", run_dir=run_dir, timeout=timeout, backend="fastp")


def trim_cutadapt(
    *,
    read1: Path,
    output1: Path,
    read2: Path | None,
    output2: Path | None,
    unpaired1: Path | None,
    unpaired2: Path | None,
    html: Path | None,
    json_report: Path | None,
    adapter: str | None,
    front: str | None,
    anywhere: str | None,
    adapter2: str | None,
    front2: str | None,
    anywhere2: str | None,
    quality_cutoff: str | None,
    minimum_length: str | None,
    maximum_length: str | None,
    max_n: str | None,
    discard_untrimmed: bool,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if unpaired1 is not None or unpaired2 is not None:
        raise MicrobiomeSuiteError(
            "--unpaired1 and --unpaired2 are only supported by --backend trimmomatic."
        )
    if html is not None:
        raise MicrobiomeSuiteError(
            "--html is only supported by --backend fastp. Use --json-report for Cutadapt."
        )
    if read2 is not None and output2 is None:
        raise MicrobiomeSuiteError("--output2 is required when --read2 is supplied.")
    if output2 is not None and read2 is None:
        raise MicrobiomeSuiteError("--read2 is required when --output2 is supplied.")
    if read2 is None and any(value is not None for value in (adapter2, front2, anywhere2)):
        raise MicrobiomeSuiteError("R2 adapter options require --read2 and --output2.")
    has_adapter = any(
        value is not None for value in (adapter, front, anywhere, adapter2, front2, anywhere2)
    )
    if (
        not any(
            value is not None
            for value in (
                adapter,
                front,
                anywhere,
                adapter2,
                front2,
                anywhere2,
                quality_cutoff,
                minimum_length,
                maximum_length,
                max_n,
            )
        )
        and not discard_untrimmed
    ):
        raise MicrobiomeSuiteError(
            "cutadapt requires at least one adapter or filtering option, such as "
            "--adapter, --front, --quality-cutoff, --minimum-length, or --max-n."
        )
    if discard_untrimmed and not has_adapter:
        raise MicrobiomeSuiteError(
            "--discard-untrimmed requires at least one adapter option, such as "
            "--adapter, --front, or --anywhere."
        )

    cutadapt = shutil.which("cutadapt")
    if cutadapt is None:
        raise MicrobiomeSuiteError(
            "Cutadapt trimming requires the external 'cutadapt' command. "
            "Install cutadapt and rerun this command."
        )

    ensure_input(read1)
    if read2 is not None:
        ensure_input(read2)
    _prepare_outputs(output1, output2, json_report, force=force)

    command = [cutadapt]
    _add_optional(command, "-a", adapter)
    _add_optional(command, "-g", front)
    _add_optional(command, "-b", anywhere)
    _add_optional(command, "-A", adapter2)
    _add_optional(command, "-G", front2)
    _add_optional(command, "-B", anywhere2)
    _add_optional(command, "-q", quality_cutoff)
    _add_optional(command, "-m", minimum_length)
    _add_optional(command, "-M", maximum_length)
    _add_optional(command, "--max-n", max_n)
    if discard_untrimmed:
        command.append("--discard-untrimmed")
    if json_report is not None:
        command.extend(["--json", str(json_report)])
    command.extend(["-j", str(threads), "-o", str(output1)])
    if read2 is not None and output2 is not None:
        command.extend(["-p", str(output2)])
    command.append(str(read1))
    if read2 is not None:
        command.append(str(read2))

    _run(command, "Cutadapt trimming failed.", run_dir=run_dir, timeout=timeout, backend="cutadapt")


def _add_optional(command: list[str], option: str, value: str | None) -> None:
    if value is not None:
        command.extend([option, value])


def _reject_options(backend: str, options: dict[str, object | None]) -> None:
    rejected = [
        option
        for option, value in options.items()
        if value is not None and value is not False and value != []
    ]
    if rejected:
        raise MicrobiomeSuiteError(f"{', '.join(rejected)} not supported by --backend {backend}.")


def trim_trimmomatic(
    *,
    read1: Path,
    output1: Path,
    read2: Path | None,
    output2: Path | None,
    unpaired1: Path | None,
    unpaired2: Path | None,
    steps: list[str],
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if not steps:
        raise MicrobiomeSuiteError(
            "trimmomatic requires at least one trimming step, such as SLIDINGWINDOW:4:20."
        )
    if read2 is not None and output2 is None:
        raise MicrobiomeSuiteError("--output2 is required when --read2 is supplied.")
    if output2 is not None and read2 is None:
        raise MicrobiomeSuiteError("--read2 is required when --output2 is supplied.")
    if read2 is not None and (unpaired1 is None or unpaired2 is None):
        raise MicrobiomeSuiteError(
            "--unpaired1 and --unpaired2 are required for paired-end Trimmomatic."
        )
    if read2 is None and (unpaired1 is not None or unpaired2 is not None):
        raise MicrobiomeSuiteError("--unpaired outputs require --read2 and --output2.")

    trimmomatic = shutil.which("trimmomatic")
    if trimmomatic is None:
        raise MicrobiomeSuiteError(
            "Trimmomatic trimming requires the external 'trimmomatic' command. "
            "Install Trimmomatic and rerun this command."
        )
    ensure_input(read1)
    if read2 is not None:
        ensure_input(read2)
    _prepare_outputs(output1, output2, unpaired1, unpaired2, force=force)

    command = [trimmomatic, "PE" if read2 is not None else "SE", "-threads", str(threads)]
    if read2 is None:
        command.extend([str(read1), str(output1)])
    else:
        command.extend(
            [
                str(read1),
                str(read2),
                str(output1),
                str(unpaired1),
                str(output2),
                str(unpaired2),
            ]
        )
    command.extend(steps)
    _run(
        command,
        "Trimmomatic trimming failed.",
        run_dir=run_dir,
        timeout=timeout,
        backend="trimmomatic",
    )


def trim_galore(
    *,
    read1: Path,
    output1: Path,
    read2: Path | None,
    output2: Path | None,
    adapter: str | None,
    adapter2: str | None,
    quality_cutoff: str | None,
    minimum_length: str | None,
    basename: str | None,
    trim_galore_version: str,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if trim_galore_version not in {"auto", "legacy", "v2"}:
        raise MicrobiomeSuiteError("--trim-galore-version must be one of: auto, legacy, v2.")
    if read2 is not None and output2 is None:
        raise MicrobiomeSuiteError("--output2 is required when --read2 is supplied.")
    if output2 is not None and read2 is None:
        raise MicrobiomeSuiteError("--read2 is required when --output2 is supplied.")
    paired = read2 is not None
    expected_output1 = _trim_galore_expected_output(output1.parent, read1, basename, paired=paired)
    if output1 != expected_output1:
        raise MicrobiomeSuiteError(
            f"Trim Galore expected output1 path is {expected_output1}; got {output1}."
        )
    if read2 is not None and output2 is not None:
        expected_output2 = _trim_galore_expected_output(
            output1.parent, read2, basename, paired=True
        )
        if output2 != expected_output2:
            raise MicrobiomeSuiteError(
                f"Trim Galore expected output2 path is {expected_output2}; got {output2}."
            )

    trim_galore_cmd = shutil.which("trim_galore")
    if trim_galore_cmd is None:
        raise MicrobiomeSuiteError(
            "Trim Galore trimming requires the external 'trim_galore' command. "
            "Install Trim Galore and rerun this command."
        )
    ensure_input(read1)
    if read2 is not None:
        ensure_input(read2)
    output_dir = output1.parent
    _prepare_outputs(output1, output2, force=force)

    command = [trim_galore_cmd]
    if read2 is not None:
        command.append("--paired")
    _add_optional(command, "--adapter", adapter)
    _add_optional(command, "--adapter2", adapter2)
    _add_optional(command, "--quality", quality_cutoff)
    _add_optional(command, "--length", minimum_length)
    command.extend(["--cores", str(threads), "--output_dir", str(output_dir)])
    _add_optional(command, "--basename", basename)
    if trim_galore_version == "v2":
        command.extend(["--engine", "v2"])
    command.append(str(read1))
    if read2 is not None:
        command.append(str(read2))
    _run(
        command,
        "Trim Galore trimming failed.",
        run_dir=run_dir,
        timeout=timeout,
        backend="trim-galore",
    )


def trim_qiime2_cutadapt(
    *,
    demux: Path,
    output: Path,
    adapter: str | None,
    front: str | None,
    anywhere: str | None,
    adapter2: str | None,
    front2: str | None,
    anywhere2: str | None,
    quality_cutoff: str | None,
    minimum_length: str | None,
    maximum_length: str | None,
    max_n: str | None,
    discard_untrimmed: bool,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    has_adapter = any(
        value is not None for value in (adapter, front, anywhere, adapter2, front2, anywhere2)
    )
    if (
        not any(
            value is not None
            for value in (
                adapter,
                front,
                anywhere,
                adapter2,
                front2,
                anywhere2,
                quality_cutoff,
                minimum_length,
                maximum_length,
                max_n,
            )
        )
        and not discard_untrimmed
    ):
        raise MicrobiomeSuiteError(
            "qiime2-cutadapt requires at least one adapter or filtering option, such as "
            "--adapter, --front, --quality-cutoff, --minimum-length, or --max-n."
        )
    if discard_untrimmed and not has_adapter:
        raise MicrobiomeSuiteError(
            "--discard-untrimmed requires at least one adapter option, such as "
            "--adapter, --front, or --anywhere."
        )

    qiime = require_qiime("QIIME 2 Cutadapt trimming")
    ensure_input(demux)
    prepare_output(output, force=force)

    paired = any(value is not None for value in (adapter2, front2, anywhere2))
    action = "trim-paired" if paired else "trim-single"
    command = [
        qiime,
        "cutadapt",
        action,
        "--i-demultiplexed-sequences",
        str(demux),
        "--o-trimmed-sequences",
        str(output),
        "--p-cores",
        str(threads),
    ]
    if paired:
        _add_optional(command, "--p-adapter-f", adapter)
        _add_optional(command, "--p-front-f", front)
        _add_optional(command, "--p-anywhere-f", anywhere)
        _add_optional(command, "--p-adapter-r", adapter2)
        _add_optional(command, "--p-front-r", front2)
        _add_optional(command, "--p-anywhere-r", anywhere2)
    else:
        _add_optional(command, "--p-adapter", adapter)
        _add_optional(command, "--p-front", front)
        _add_optional(command, "--p-anywhere", anywhere)
    _add_qiime_quality_cutoff(command, quality_cutoff)
    _add_optional(command, "--p-minimum-length", minimum_length)
    _add_optional(command, "--p-maximum-length", maximum_length)
    _add_optional(command, "--p-max-n", max_n)
    if discard_untrimmed:
        command.append("--p-discard-untrimmed")

    run_qiime(
        command,
        "QIIME 2 Cutadapt trimming failed.",
        run_dir=run_dir,
        timeout=timeout,
        task="trim",
        backend="qiime2-cutadapt",
    )


def _add_qiime_quality_cutoff(command: list[str], quality_cutoff: str | None) -> None:
    if quality_cutoff is None:
        return
    values = quality_cutoff.split(",", maxsplit=1)
    if len(values) == 1:
        command.extend(["--p-quality-cutoff-3end", values[0]])
        return
    command.extend(["--p-quality-cutoff-5end", values[0]])
    command.extend(["--p-quality-cutoff-3end", values[1]])


def _trim_galore_expected_output(
    output_dir: Path, read: Path, basename: str | None, *, paired: bool
) -> Path:
    if basename is not None:
        if paired:
            suffix = "_val_2.fq.gz" if "_R2" in read.name or "_2" in read.name else "_val_1.fq.gz"
            return output_dir / f"{basename}{suffix}"
        return output_dir / f"{basename}.fq.gz"
    stem = read.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    if paired:
        suffix = "_val_2.fq.gz" if "_R2" in read.name or "_2" in read.name else "_val_1.fq.gz"
    else:
        suffix = "_trimmed.fq.gz"
    return output_dir / f"{stem}{suffix}"


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
        log=CommandLog(task="trim", backend=backend),
    )


def _prepare_outputs(*outputs: Path | None, force: bool) -> None:
    for output in outputs:
        if output is not None:
            prepare_output(output, force=force)
