"""End-to-end mothur MiSeq SOP: FASTQ directory to OTU table and taxonomy."""

from __future__ import annotations

import re
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import prepare_output
from microsuite.methods.cluster import cluster
from microsuite.methods.mothur import run_mothur, select_output
from microsuite.methods.tax_classify import tax_classify

# Matches common FASTQ mate-pair naming schemes, case-insensitively:
#   sampleA_R1.fastq.gz
#   sampleA_S1_L001_R1_001.fastq.gz   (bcl2fastq default output; _S<n>/_L<lane>
#                                      are stripped from the sample name)
#   sampleA_r1.fastq.gz               (lowercase mate marker)
#   sampleA.R1.fastq.gz               (dot separator)
_MATE = re.compile(
    r"^(?P<sample>.+?)(?:_S\d+)?(?:_L\d{3})?[._](?:R)?(?P<mate>[12])(?:_001)?"
    r"\.f(?:ast)?q(?:\.gz)?$",
    re.IGNORECASE,
)
# Any file that "looks like" a FASTQ file, regardless of whether _MATE can parse it.
_FASTQ_SUFFIX = re.compile(r"\.f(?:ast)?q(?:\.gz)?$", re.IGNORECASE)


def write_stability_file(reads_dir: Path, output: Path) -> Path:
    """Write mothur's stability file: sample name, R1 path, R2 path per line."""
    pairs: dict[str, dict[str, Path]] = {}
    for path in sorted(reads_dir.iterdir()):
        match = _MATE.match(path.name)
        if match is None:
            if _FASTQ_SUFFIX.search(path.name):
                raise MicrobiomeSuiteError(
                    f"Cannot determine sample/mate for FASTQ file '{path.name}'. "
                    "Expected a name like 'sampleA_R1.fastq.gz' (or "
                    "'sampleA_S1_L001_R1_001.fastq.gz' for Illumina bcl2fastq output)."
                )
            continue
        sample = match.group("sample")
        mate = match.group("mate")
        existing = pairs.setdefault(sample, {}).get(mate)
        if existing is not None:
            raise MicrobiomeSuiteError(
                f"Multiple files map to sample '{sample}' mate {mate}: "
                f"'{existing.name}' and '{path.name}'. mothur's stability file "
                "needs exactly one R1 and one R2 per sample; this usually means "
                "reads from more than one sequencing lane were not concatenated "
                "first. Concatenate each sample's lanes into a single R1/R2 pair "
                "before running this workflow."
            )
        pairs[sample][mate] = path

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

    # mothur derives output names from the stability file's stem ("stability"),
    # so the assembled-contigs path is deterministic before make.contigs runs.
    # Guarding it here -- like cluster() and tax_classify() guard their own
    # deliverables -- makes a rerun without --force fail before overwriting the
    # prior contigs, instead of after.
    contigs_fasta = work_dir / "stability.trim.contigs.fasta"
    prepare_output(contigs_fasta, force=force)

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
    # make.contigs's count table is the ONLY carrier of read->sample identity --
    # mothur does not rename reads to encode a sample. Without it, cluster()
    # falls back to a group-less dereplication and every downstream table ends
    # up with a single sample column regardless of how many samples were fed in.
    contigs_count_table = select_output(outputs, ".count_table", step="make.contigs")

    output_table = output_dir / "table.tsv"
    output_rep_seqs = output_dir / "rep-seqs.fasta"
    output_otu_list = output_dir / "otu.list"
    output_count_table = output_dir / "table.count_table"
    cluster(
        backend="mothur",
        rep_seqs=contigs,
        output_table=output_table,
        output_rep_seqs=output_rep_seqs,
        reference_alignment=reference_alignment,
        count_table=contigs_count_table,
        identity=identity,
        output_otu_list=output_otu_list,
        output_count_table=output_count_table,
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
        otu_list=output_otu_list,
        count_table=output_count_table,
        force=force,
        run_dir=run_dir / "taxonomy",
        timeout=timeout,
    )
