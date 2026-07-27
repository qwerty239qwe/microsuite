"""End-to-end mothur MiSeq SOP: FASTQ directory to OTU table and taxonomy."""

from __future__ import annotations

import re
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.cluster import cluster
from microsuite.methods.mothur import run_mothur, select_output
from microsuite.methods.tax_classify import tax_classify

_MATE = re.compile(r"^(?P<sample>.+?)_(?:R)?(?P<mate>[12])(?:_001)?\.f(?:ast)?q(?:\.gz)?$")


def write_stability_file(reads_dir: Path, output: Path) -> Path:
    """Write mothur's stability file: sample name, R1 path, R2 path per line."""
    pairs: dict[str, dict[str, Path]] = {}
    for path in sorted(reads_dir.iterdir()):
        match = _MATE.match(path.name)
        if match is None:
            continue
        pairs.setdefault(match.group("sample"), {})[match.group("mate")] = path

    if not pairs:
        raise MicrobiomeSuiteError(f"No paired FASTQ files found in {reads_dir}.")

    unpaired = sorted(sample for sample, mates in pairs.items() if len(mates) != 2)
    if unpaired:
        raise MicrobiomeSuiteError(
            f"Samples missing a mate file: {', '.join(unpaired)}. "
            "mothur's make.contigs requires both reads of every pair."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        for sample in sorted(pairs):
            mates = pairs[sample]
            handle.write(f"{sample}\t{mates['1']}\t{mates['2']}\n")
    return output


def run_mothur_sop(
    *,
    reads_dir: Path,
    output_dir: Path,
    reference_alignment: Path,
    taxonomy_reference: Path,
    taxonomy_map: Path,
    identity: float = 0.97,
    force: bool = False,
    timeout: float | None = None,
) -> None:
    """Run make.contigs, then the mothur cluster and taxonomy backends."""
    if not reads_dir.is_dir():
        raise MicrobiomeSuiteError(f"Reads directory does not exist: {reads_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    run_dir = output_dir / "logs"

    work_dir = output_dir / "contigs"
    stability = write_stability_file(reads_dir, work_dir / "stability.files")

    outputs = run_mothur(
        "make.contigs",
        {"file": str(stability)},
        work_dir=work_dir,
        run_dir=run_dir / "make_contigs",
        timeout=timeout,
    )
    # make.contigs emits both .trim.contigs.fasta and .scrap.contigs.fasta;
    # the scrap file holds the reads that failed assembly.
    contigs = select_output(outputs, ".fasta", step="make.contigs", exclude=("scrap",))

    output_table = output_dir / "table.tsv"
    output_rep_seqs = output_dir / "rep-seqs.fasta"
    cluster(
        backend="mothur",
        rep_seqs=contigs,
        output_table=output_table,
        output_rep_seqs=output_rep_seqs,
        reference_alignment=reference_alignment,
        identity=identity,
        force=force,
        run_dir=run_dir / "cluster",
        timeout=timeout,
    )

    tax_classify(
        backend="mothur",
        rep_seqs=output_rep_seqs,
        output=output_dir / "taxonomy.tsv",
        taxonomy_reference=taxonomy_reference,
        taxonomy_map=taxonomy_map,
        force=force,
        run_dir=run_dir / "taxonomy",
        timeout=timeout,
    )
