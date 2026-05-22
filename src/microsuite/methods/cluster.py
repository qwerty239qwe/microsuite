from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.runtime.runner import CommandLog, run_command

SUPPORTED_BACKENDS = ("vsearch", "usearch")


def cluster(
    *,
    backend: str,
    rep_seqs: Path,
    output_table: Path,
    output_rep_seqs: Path,
    table: Path | None = None,
    identity: float = 0.97,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    if backend == "vsearch":
        if table is None:
            raise MicrobiomeSuiteError("--table is required for --backend vsearch.")
        cluster_vsearch(
            table=table,
            rep_seqs=rep_seqs,
            output_table=output_table,
            output_rep_seqs=output_rep_seqs,
            identity=identity,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "usearch":
        cluster_usearch(
            rep_seqs=rep_seqs,
            output_uc=output_table,
            output_centroids=output_rep_seqs,
            identity=identity,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    else:
        backends = ", ".join(SUPPORTED_BACKENDS)
        raise MicrobiomeSuiteError(
            f"Unsupported cluster backend '{backend}'. Choose one of: {backends}"
        )


def cluster_vsearch(
    *,
    table: Path,
    rep_seqs: Path,
    output_table: Path,
    output_rep_seqs: Path,
    identity: float,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
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
    run_command(
        command,
        "QIIME 2 VSEARCH clustering failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="cluster", backend="vsearch"),
    )


def cluster_usearch(
    *,
    rep_seqs: Path,
    output_uc: Path,
    output_centroids: Path,
    identity: float,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if not 0 < identity <= 1:
        raise MicrobiomeSuiteError("--identity must be greater than 0 and less than or equal to 1.")
    usearch = shutil.which("usearch") or shutil.which("usearch12")
    if usearch is None:
        raise MicrobiomeSuiteError(
            "USEARCH clustering requires the external 'usearch' command. "
            "Install USEARCH 12 or use the microsuite/usearch container and rerun this command."
        )

    ensure_input(rep_seqs)
    prepare_output(output_uc, force=force)
    prepare_output(output_centroids, force=force)

    command = [
        usearch,
        "-cluster_fast",
        str(rep_seqs),
        "-id",
        str(identity),
        "-centroids",
        str(output_centroids),
        "-uc",
        str(output_uc),
    ]
    run_command(
        command,
        "USEARCH clustering failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="cluster", backend="usearch"),
    )
