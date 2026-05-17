from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output

SUPPORTED_METHODS = ("qiime2", "kraken2", "dada2")
PLANNED_METHODS = ("kraken2", "dada2")


def tax_classify(
    *,
    backend: str,
    rep_seqs: Path,
    output: Path,
    classifier: Path | None = None,
    threads: int = 1,
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend == "qiime2":
        tax_classify_qiime2(
            rep_seqs=rep_seqs,
            classifier=classifier,
            output=output,
            threads=threads,
            force=force,
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


def tax_classify_qiime2(
    *,
    rep_seqs: Path,
    classifier: Path | None,
    output: Path,
    threads: int,
    force: bool,
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
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip() or "QIIME 2 classification failed."
        raise MicrobiomeSuiteError(message)
