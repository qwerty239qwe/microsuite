"""End-to-end smoke test for the mothur MiSeq SOP workflow.

Gated behind ``MICROSUITE_RUN_MOTHUR_SMOKE=1`` and a mothur binary on PATH,
because it runs the real 12-command pipeline.

**Why this test exists.** The mothur backends shipped with a fully green unit
suite that mocked every subprocess call, and five defects survived it — all of
them producing a complete, well-formed, *wrong* result rather than an error:

* the sample→read mapping was dropped, so every table had one column
* ``screen.seqs`` was given the wrong output extension and its count table is
  emitted only when sequences are actually removed
* ``chimera.vsearch``'s outputs invert depending on whether the count table
  carries sample groups
* per-OTU taxonomy was consensused over a fraction of each OTU's members

None of those are visible when subprocess is mocked, because the mock returns
whatever stdout the test author *believed* mothur produces. This test runs the
real thing and asserts on the biology.

The dataset is generated deterministically rather than committed, so there is
no binary fixture to drift.
"""

from __future__ import annotations

import os
import random
import shutil
from pathlib import Path

import pytest

from microsuite.workflows.mothur_sop import run_mothur_sop

pytestmark = [
    pytest.mark.skipif(
        os.environ.get("MICROSUITE_RUN_MOTHUR_SMOKE") != "1",
        reason="set MICROSUITE_RUN_MOTHUR_SMOKE=1 to run the real mothur pipeline",
    ),
    pytest.mark.skipif(shutil.which("mothur") is None, reason="mothur is not installed"),
]

# Read counts per sample per species. Deliberately asymmetric: a bug that
# collapses samples together, or that loses one, changes these numbers.
DESIGN = {
    "sampleA": {"spA": 30, "spB": 12, "spC": 3},
    "sampleB": {"spA": 6, "spB": 28, "spC": 20},
}
_RC = str.maketrans("ACGT", "TGCA")


def _build_dataset(root: Path) -> None:
    """Write paired FASTQ for two samples, an aligned reference, and a trainset."""
    rng = random.Random(42)
    core = "".join(rng.choice("ACGT") for _ in range(250))

    def mutate(seq: str, count: int) -> str:
        chars = list(seq)
        for position in rng.sample(range(30, 220), count):
            chars[position] = rng.choice([c for c in "ACGT" if c != chars[position]])
        return "".join(chars)

    species = {"spA": core, "spB": mutate(core, 25), "spC": mutate(core, 40)}

    reads = root / "reads"
    reads.mkdir(parents=True, exist_ok=True)
    quality = "I" * 150
    for sample, composition in DESIGN.items():
        r1: list[str] = []
        r2: list[str] = []
        index = 0
        for name, count in composition.items():
            for _ in range(count):
                chars = list(species[name])
                # A little sequencing error, so pre.cluster and chimera.vsearch
                # have something to do rather than seeing identical reads.
                for position in rng.sample(range(len(chars)), 2):
                    chars[position] = rng.choice("ACGT")
                seq = "".join(chars)
                read_id = f"{sample}:{index}"
                index += 1
                # 150bp mates over a 250bp amplicon leaves a 50bp overlap for
                # make.contigs to assemble on.
                r1.append(f"@{read_id} 1:N:0\n{seq[:150]}\n+\n{quality}\n")
                r2.append(f"@{read_id} 2:N:0\n{seq[100:].translate(_RC)[::-1]}\n+\n{quality}\n")
        (reads / f"{sample}_R1.fastq").write_text("".join(r1), encoding="utf-8")
        (reads / f"{sample}_R2.fastq").write_text("".join(r2), encoding="utf-8")

    # align.seqs rejects an unaligned reference, so the gap column is required.
    (root / "ref.align").write_text(
        f">refA\n{core[:120] + '-' * 15 + core[120:]}\n", encoding="utf-8"
    )
    (root / "train.fasta").write_text(
        "".join(f">{k}\n{v}\n" for k, v in species.items()), encoding="utf-8"
    )
    (root / "train.tax").write_text(
        "spA\tBacteria;Firmicutes;Bacilli;\n"
        "spB\tBacteria;Bacteroidota;Bacteroidia;\n"
        "spC\tBacteria;Proteobacteria;Gammaproteobacteria;\n",
        encoding="utf-8",
    )


def _read_tsv(path: Path) -> list[list[str]]:
    return [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines() if line]


@pytest.fixture(scope="module")
def sop_output(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Run the full workflow once; every assertion below reads its outputs."""
    root = tmp_path_factory.mktemp("mothur-smoke")
    _build_dataset(root)
    output_dir = root / "out"
    run_mothur_sop(
        reads_dir=root / "reads",
        output_dir=output_dir,
        reference_alignment=root / "ref.align",
        taxonomy_reference=root / "train.fasta",
        taxonomy_map=root / "train.tax",
    )
    return output_dir


def test_workflow_produces_every_declared_output(sop_output: Path) -> None:
    for name in ("table.tsv", "table.shared", "rep-seqs.fasta", "taxonomy.tsv"):
        path = sop_output / name
        assert path.exists(), f"{name} was not written"
        assert path.stat().st_size > 0, f"{name} is empty"


def test_table_keeps_every_sample_as_its_own_column(sop_output: Path) -> None:
    # The regression that motivated this file: make.contigs' count table is the
    # only carrier of read->sample identity, and dropping it yielded a
    # single-column table with no error raised anywhere.
    header, *rows = _read_tsv(sop_output / "table.tsv")

    assert header[0] == "feature-id"
    assert header[1:] == sorted(DESIGN), f"expected one column per sample, got {header[1:]}"
    assert rows, "table has no feature rows"


def test_table_recovers_the_designed_community_structure(sop_output: Path) -> None:
    header, *rows = _read_tsv(sop_output / "table.tsv")
    samples = header[1:]
    per_sample_total = {
        sample: sum(int(row[i + 1]) for row in rows) for i, sample in enumerate(samples)
    }

    # Clustering at 97% should recover one OTU per species, not one per read.
    assert len(rows) == 3, f"expected 3 OTUs for 3 species, got {len(rows)}"

    # Every read should survive to the table. Allow a small loss to chimera and
    # quality filtering, but not the ~50% a dropped sample or lane would cost.
    for sample, expected in DESIGN.items():
        assert per_sample_total[sample] >= sum(expected.values()) * 0.8, (
            f"{sample} kept {per_sample_total[sample]} of {sum(expected.values())} reads"
        )

    # The two samples have deliberately different compositions; identical
    # columns would mean sample identity collapsed somewhere upstream.
    columns = [tuple(row[i + 1] for row in rows) for i in range(len(samples))]
    assert columns[0] != columns[1], "sample columns are identical; grouping was lost"


def test_taxonomy_is_per_otu_consensus_over_full_membership(sop_output: Path) -> None:
    # classify.otu must see every sequence in each OTU, not just its
    # representative. When it saw only the representatives, the Size column
    # reported 29/2/1 against true OTU totals of 40/34/23 -- a plausible,
    # well-formed, wrong consensus.
    shared_header, *shared_rows = _read_tsv(sop_output / "table.shared")
    otus = shared_header[3:]
    shared_total = {otu: sum(int(row[3 + i]) for row in shared_rows) for i, otu in enumerate(otus)}

    tax_header, *tax_rows = _read_tsv(sop_output / "taxonomy.tsv")
    assert tax_header[:3] == ["OTU", "Size", "Taxonomy"], (
        f"expected classify.otu consensus output, got {tax_header}"
    )

    for otu, size, _lineage in (row[:3] for row in tax_rows):
        assert int(size) == shared_total[otu], (
            f"{otu}: consensus covered {size} sequences but the OTU holds "
            f"{shared_total[otu]} -- classification saw only part of the OTU"
        )


def test_clustering_actually_clustered(sop_output: Path) -> None:
    # A degenerate distance matrix makes mothur emit a .list where every unique
    # sequence is its own OTU. Every later step then succeeds, so the run exits
    # 0 with a full table and full taxonomy that happen to describe no
    # clustering at all. mothur exits 1 on the two ways we know to provoke this,
    # so run_command catches it -- this asserts the outcome directly rather than
    # relying on that staying true.
    _header, *rows = _read_tsv(sop_output / "table.tsv")
    unique_seqs = sum(
        1
        for line in (sop_output / "unique-seqs.fasta").read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    )

    assert len(rows) < unique_seqs, (
        f"{len(rows)} OTUs from {unique_seqs} unique sequences -- clustering "
        "collapsed nothing, which is what a blank distance matrix produces"
    )


def test_representative_sequences_are_unaligned(sop_output: Path) -> None:
    # get.oturep reads the filtered alignment, so its output carries gap
    # characters unless degapped. The vsearch and usearch backends return
    # unaligned centroids under this same option, and downstream consumers
    # (mafft-fasttree, qiime2 classifiers) assume unaligned input.
    sequences = [
        line
        for line in (sop_output / "rep-seqs.fasta").read_text(encoding="utf-8").splitlines()
        if line and not line.startswith(">")
    ]

    assert sequences, "rep-seqs.fasta has no sequence lines"
    for line in sequences:
        assert "-" not in line and "." not in line, (
            f"rep-seqs.fasta still carries alignment gaps: {line[:60]}"
        )


def test_taxonomy_assigns_a_real_lineage_to_every_otu(sop_output: Path) -> None:
    # Mismatched name spaces between classify.seqs and classify.otu produce a
    # complete table of "unknown" that looks entirely normal.
    _header, *rows = _read_tsv(sop_output / "taxonomy.tsv")

    assert rows, "taxonomy table has no rows"
    for otu, _size, lineage in (row[:3] for row in rows):
        assert "Bacteria" in lineage, f"{otu} was not classified: {lineage}"
        assert "unknown" not in lineage.lower(), f"{otu} came back unclassified: {lineage}"
