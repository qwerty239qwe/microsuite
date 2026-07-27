from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.runner import CommandLog, run_command

SUPPORTED_BACKENDS = ("vsearch", "usearch", "qiime2-vsearch")


def cluster(
    *,
    backend: str,
    rep_seqs: Path,
    output_table: Path,
    output_rep_seqs: Path,
    table: Path | None = None,
    reads: Path | None = None,
    identity: float = 0.97,
    output_uc: Path | None = None,
    sample_delimiter: str = "_",
    sample_field: int = 0,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "cluster")
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
            reads=reads,
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
    reads: Path | None = None,
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
    if reads is not None:
        ensure_input(reads)
    prepare_output(output_table, force=force)
    prepare_output(output_centroids, force=force)
    output_uc = output_uc or output_table.with_suffix(".uc")
    prepare_output(output_uc, force=force)

    run_command(
        [
            vsearch,
            "--cluster_fast",
            str(rep_seqs),
            "--id",
            str(identity),
            "--centroids",
            str(output_centroids),
            "--uc",
            str(output_uc),
        ],
        "VSEARCH clustering failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(task="cluster", backend="vsearch"),
    )

    # The clustering .uc only describes the (possibly dereplicated) rep_seqs, so
    # counting it would mis-attribute per-sample abundances. When the full read
    # set is supplied, map every read back to the centroids and build the table
    # from that mapping instead, matching the standard derep -> cluster -> remap
    # OTU workflow.
    if reads is not None:
        map_uc = output_table.with_suffix(".map.uc")
        prepare_output(map_uc, force=force)
        run_command(
            [
                vsearch,
                "--usearch_global",
                str(reads),
                "--db",
                str(output_centroids),
                "--id",
                str(identity),
                "--uc",
                str(map_uc),
            ],
            "VSEARCH read mapping failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(task="cluster", backend="vsearch"),
        )
        table_uc = map_uc
    else:
        table_uc = output_uc

    write_otu_table_from_uc(
        table_uc,
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
            size = _size_from_label(fields[8])
            samples.add(sample)
            counts.setdefault(otu, {})
            counts[otu][sample] = counts[otu].get(sample, 0) + size

    ordered_samples = sorted(samples)
    with output.open("w", encoding="utf-8", newline="") as handle:
        handle.write("feature-id\t" + "\t".join(ordered_samples) + "\n")
        for otu in sorted(counts):
            values = [str(counts[otu].get(sample, 0)) for sample in ordered_samples]
            handle.write(otu + "\t" + "\t".join(values) + "\n")


def _clean_uc_label(label: str) -> str:
    return label.split(";", 1)[0]


def _size_from_label(label: str) -> int:
    # vsearch annotates dereplicated abundance as ';size=N'; absent it, each
    # .uc record represents a single read.
    for part in label.split(";"):
        if part.startswith("size="):
            try:
                return int(part[len("size=") :])
            except ValueError:
                return 1
    return 1


def _sample_from_label(label: str, *, sample_delimiter: str, sample_field: int) -> str:
    if sample_field < 0:
        raise MicrobiomeSuiteError("--sample-field must be zero or greater.")
    if not sample_delimiter:
        return label
    fields = label.split(sample_delimiter)
    if sample_field >= len(fields):
        return label
    return fields[sample_field]


def write_otu_table_from_shared(shared: Path, output: Path) -> None:
    """Convert mothur's sample-major .shared into microsuite's feature-major TSV.

    .shared columns are: label, Group, numOtus, then one column per OTU.
    """
    rows = [
        line.rstrip("\n").split("\t")
        for line in shared.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if len(rows) < 2:
        raise MicrobiomeSuiteError(f"mothur .shared file has no rows: {shared}")

    header, *records = rows
    otus = header[3:]
    samples = [record[1] for record in records]

    with output.open("w", encoding="utf-8", newline="") as handle:
        handle.write("feature-id\t" + "\t".join(samples) + "\n")
        for index, otu in enumerate(otus):
            counts = [record[3 + index] for record in records]
            handle.write(otu + "\t" + "\t".join(counts) + "\n")
