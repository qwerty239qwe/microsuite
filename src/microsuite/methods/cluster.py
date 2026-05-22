from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.runtime.runner import CommandLog, run_command

SUPPORTED_BACKENDS = ("vsearch", "usearch", "qiime2-vsearch")


def cluster(
    *,
    backend: str,
    rep_seqs: Path,
    output_table: Path,
    output_rep_seqs: Path,
    table: Path | None = None,
    identity: float = 0.97,
    output_uc: Path | None = None,
    sample_delimiter: str = "_",
    sample_field: int = 0,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    if backend == "qiime2-vsearch":
        if table is None:
            raise MicrobiomeSuiteError("--table is required for --backend qiime2-vsearch.")
        cluster_qiime2_vsearch(
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
    if backend == "vsearch":
        cluster_vsearch(
            rep_seqs=rep_seqs,
            output_table=output_table,
            output_centroids=output_rep_seqs,
            output_uc=output_uc,
            identity=identity,
            sample_delimiter=sample_delimiter,
            sample_field=sample_field,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "usearch":
        cluster_usearch(
            rep_seqs=rep_seqs,
            output_table=output_table,
            output_centroids=output_rep_seqs,
            output_uc=output_uc,
            identity=identity,
            sample_delimiter=sample_delimiter,
            sample_field=sample_field,
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


def cluster_qiime2_vsearch(
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
        log=CommandLog(task="cluster", backend="qiime2-vsearch"),
    )


def cluster_vsearch(
    *,
    rep_seqs: Path,
    output_table: Path,
    output_centroids: Path,
    output_uc: Path | None,
    identity: float,
    sample_delimiter: str,
    sample_field: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if not 0 < identity <= 1:
        raise MicrobiomeSuiteError("--identity must be greater than 0 and less than or equal to 1.")
    vsearch = shutil.which("vsearch")
    if vsearch is None:
        raise MicrobiomeSuiteError(
            "VSEARCH clustering requires the external 'vsearch' command. "
            "Install VSEARCH or use the microsuite/vsearch container and rerun this command."
        )

    ensure_input(rep_seqs)
    prepare_output(output_table, force=force)
    prepare_output(output_centroids, force=force)
    output_uc = output_uc or output_table.with_suffix(".uc")
    prepare_output(output_uc, force=force)

    command = [
        vsearch,
        "--cluster_fast",
        str(rep_seqs),
        "--id",
        str(identity),
        "--centroids",
        str(output_centroids),
        "--uc",
        str(output_uc),
    ]
    run_command(
        command,
        "VSEARCH clustering failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="cluster", backend="vsearch"),
    )
    write_otu_table_from_uc(
        output_uc,
        output_table,
        sample_delimiter=sample_delimiter,
        sample_field=sample_field,
    )


def cluster_usearch(
    *,
    rep_seqs: Path,
    output_table: Path,
    output_centroids: Path,
    output_uc: Path | None,
    identity: float,
    sample_delimiter: str,
    sample_field: int,
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
    prepare_output(output_table, force=force)
    prepare_output(output_centroids, force=force)
    output_uc = output_uc or output_table.with_suffix(".uc")
    prepare_output(output_uc, force=force)

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
    write_otu_table_from_uc(
        output_uc,
        output_table,
        sample_delimiter=sample_delimiter,
        sample_field=sample_field,
    )


def write_otu_table_from_uc(
    uc: Path,
    output: Path,
    *,
    sample_delimiter: str,
    sample_field: int,
) -> None:
    counts: dict[str, dict[str, int]] = {}
    samples: set[str] = set()
    with uc.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10 or fields[0] not in {"S", "H"}:
                continue
            query = _clean_uc_label(fields[8])
            target = _clean_uc_label(fields[9])
            otu = query if target == "*" else target
            sample = _sample_from_label(
                query, sample_delimiter=sample_delimiter, sample_field=sample_field
            )
            samples.add(sample)
            counts.setdefault(otu, {})
            counts[otu][sample] = counts[otu].get(sample, 0) + 1

    ordered_samples = sorted(samples)
    with output.open("w", encoding="utf-8", newline="") as handle:
        handle.write("feature-id\t" + "\t".join(ordered_samples) + "\n")
        for otu in sorted(counts):
            values = [str(counts[otu].get(sample, 0)) for sample in ordered_samples]
            handle.write(otu + "\t" + "\t".join(values) + "\n")


def _clean_uc_label(label: str) -> str:
    return label.split(";", 1)[0]


def _sample_from_label(label: str, *, sample_delimiter: str, sample_field: int) -> str:
    if sample_field < 0:
        raise MicrobiomeSuiteError("--sample-field must be zero or greater.")
    if not sample_delimiter:
        return label
    fields = label.split(sample_delimiter)
    if sample_field >= len(fields):
        return label
    return fields[sample_field]
