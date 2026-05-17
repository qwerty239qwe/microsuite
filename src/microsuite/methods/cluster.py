from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output

SUPPORTED_BACKENDS = ("vsearch",)


def cluster(
    *,
    backend: str,
    table: Path,
    rep_seqs: Path,
    output_table: Path,
    output_rep_seqs: Path,
    identity: float = 0.97,
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend != "vsearch":
        backends = ", ".join(SUPPORTED_BACKENDS)
        raise MicrobiomeSuiteError(
            f"Unsupported cluster backend '{backend}'. Choose one of: {backends}"
        )
    cluster_vsearch(
        table=table,
        rep_seqs=rep_seqs,
        output_table=output_table,
        output_rep_seqs=output_rep_seqs,
        identity=identity,
        force=force,
    )


def cluster_vsearch(
    *,
    table: Path,
    rep_seqs: Path,
    output_table: Path,
    output_rep_seqs: Path,
    identity: float,
    force: bool,
) -> None:
    if not 0 < identity <= 1:
        raise MicrobiomeSuiteError("--identity must be greater than 0 and less than or equal to 1.")
    qiime = shutil.which("qiime")
    if qiime is None:
        raise MicrobiomeSuiteError(
            "VSEARCH clustering requires the external 'qiime' command with the vsearch plugin. "
            "Activate a QIIME 2 environment and rerun this command."
        )

    ensure_input(table)
    ensure_input(rep_seqs)
    prepare_output(output_table, force=force)
    prepare_output(output_rep_seqs, force=force)

    command = [
        qiime,
        "vsearch",
        "cluster-features-de-novo",
        "--i-table",
        str(table),
        "--i-sequences",
        str(rep_seqs),
        "--p-perc-identity",
        str(identity),
        "--o-clustered-table",
        str(output_table),
        "--o-clustered-sequences",
        str(output_rep_seqs),
    ]
    result = subprocess.run(command, check=False, text=True, capture_output=True)
    if result.returncode != 0:
        message = result.stderr.strip() or result.stdout.strip()
        message = message or "QIIME 2 VSEARCH clustering failed."
        raise MicrobiomeSuiteError(message)
