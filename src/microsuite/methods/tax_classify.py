from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_METHODS = ("qiime2", "kraken2", "dada2")
PLANNED_METHODS = ("dada2",)


def tax_classify(
    *,
    backend: str,
    rep_seqs: Path,
    output: Path,
    classifier: Path | None = None,
    threads: int | str = 1,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    resolved_threads = resolve_threads(threads)
    if backend == "qiime2":
        tax_classify_qiime2(
            rep_seqs=rep_seqs,
            classifier=classifier,
            output=output,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "kraken2":
        tax_classify_kraken2(
            reads=rep_seqs,
            database=classifier,
            output=output,
            threads=resolved_threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend in PLANNED_METHODS:
        raise MicrobiomeSuiteError(
            f"Taxonomy classification backend '{backend}' is registered but not implemented yet. "
            "Use --backend qiime2 for now."
        )
    raise MicrobiomeSuiteError(
        f"Unsupported taxonomy classification backend '{backend}'. "
        f"Choose one of: {', '.join(SUPPORTED_METHODS)}"
    )


def tax_classify_kraken2(
    *,
    reads: Path,
    database: Path | None,
    output: Path,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if database is None:
        raise MicrobiomeSuiteError("--classifier is required for --backend kraken2.")
    kraken2 = shutil.which("kraken2")
    if kraken2 is None:
        raise MicrobiomeSuiteError(
            "Kraken2 taxonomy classification requires the external 'kraken2' command. "
            "Install Kraken2 or use the microsuite/kraken2 container and rerun this command."
        )

    ensure_input(reads)
    if not database.exists():
        raise MicrobiomeSuiteError(f"Kraken2 database does not exist: {database}")
    prepare_output(output, force=force)
    per_read_output = output.with_suffix(".kraken")
    prepare_output(per_read_output, force=force)

    command = [
        kraken2,
        "--db",
        str(database),
        "--threads",
        str(threads),
        "--report",
        str(output),
        "--output",
        str(per_read_output),
        str(reads),
    ]
    run_command(
        command,
        "Kraken2 taxonomy classification failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="tax_classify",
            backend="kraken2",
            outputs={"report": str(output), "per_read": str(per_read_output)},
        ),
    )


def tax_classify_qiime2(
    *,
    rep_seqs: Path,
    classifier: Path | None,
    output: Path,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if classifier is None:
        raise MicrobiomeSuiteError("--classifier is required for --backend qiime2.")
    qiime = shutil.which("qiime")
    if qiime is None:
        raise MicrobiomeSuiteError(
            "QIIME 2 taxonomy classification requires the external 'qiime' command. "
            "Activate a QIIME 2 environment and rerun this command."
        )

    ensure_input(rep_seqs)
    ensure_input(classifier)
    prepare_output(output, force=force)

    command = [
        qiime,
        "feature-classifier",
        "classify-sklearn",
        "--i-classifier",
        str(classifier),
        "--i-reads",
        str(rep_seqs),
        "--o-classification",
        str(output),
        "--p-n-jobs",
        str(threads),
    ]
    run_command(
        command,
        "QIIME 2 classification failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="tax_classify", backend="qiime2"),
    )
