from __future__ import annotations

from pathlib import Path

import pytest

import microsuite.api as api
from microsuite._errors import MicrobiomeSuiteError
from microsuite.api import DataAsset, MicrobiomeDataset, read_table
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def test_dataset_table_roundtrip(tmp_path: Path) -> None:
    table = tmp_path / "table.h5ad"
    dataset = MicrobiomeDataset(fixture_adata())

    dataset.write(table)
    reloaded = MicrobiomeDataset.read(table)

    assert reloaded.adata is not None
    assert reloaded.adata.shape == dataset.adata.shape
    assert read_table(table).shape == dataset.adata.shape


def test_dataset_registers_path_backed_sequence_and_read_assets(tmp_path: Path) -> None:
    fasta = tmp_path / "rep-seqs.fasta"
    r1 = tmp_path / "S1_R1.fastq.gz"
    r2 = tmp_path / "S1_R2.fastq.gz"
    tree = tmp_path / "tree.nwk"
    classifier = tmp_path / "classifier.qza"
    for path in (fasta, r1, r2, tree, classifier):
        path.touch()

    dataset = MicrobiomeDataset()
    rep_seqs = dataset.add_sequences("rep_seqs", fasta)
    reads = dataset.add_reads("S1", r1=r1, r2=r2)
    rooted_tree = dataset.add_tree("rooted_tree", tree)
    silva = dataset.add_reference("silva", classifier, kind="classifier")

    assert rep_seqs == DataAsset(
        name="rep_seqs", kind="fasta", path=fasta, format=".fasta"
    )
    assert reads.kind == "paired_fastq"
    assert reads.sample_id == "S1"
    assert reads.path == r1
    assert reads.paired_path == r2
    assert reads.format == ".fastq.gz"
    assert rooted_tree.kind == "tree"
    assert silva.kind == "classifier"
    assert dataset.list_assets(kind="paired_fastq") == [reads]
    assert dataset.get_asset("rep_seqs") is rep_seqs


def test_dataset_registers_qiime_artifact_and_future_reference_dir(
    tmp_path: Path,
) -> None:
    table_qza = tmp_path / "table.qza"
    kraken_db = tmp_path / "kraken_db"
    table_qza.touch()
    kraken_db.mkdir()

    dataset = MicrobiomeDataset()
    qza = dataset.add_artifact("qiime_table", table_qza)
    reference = dataset.add_reference("kraken", kraken_db)

    assert qza.kind == "qza"
    assert qza.format == ".qza"
    assert reference.kind == "reference_db"
    assert reference.format == ""


def test_dataset_asset_validation_reports_missing_paths(tmp_path: Path) -> None:
    dataset = MicrobiomeDataset()
    dataset.add_reads("S1", r1=tmp_path / "missing.fastq.gz", require_exists=False)

    with pytest.raises(MicrobiomeSuiteError, match="does not exist"):
        dataset.validate_assets()


def test_dataset_rejects_duplicate_names_and_wrong_formats(tmp_path: Path) -> None:
    fasta = tmp_path / "rep-seqs.fasta"
    not_fasta = tmp_path / "rep-seqs.txt"
    fasta.touch()
    not_fasta.touch()
    dataset = MicrobiomeDataset()

    dataset.add_sequences("rep_seqs", fasta)

    with pytest.raises(MicrobiomeSuiteError, match="already exists"):
        dataset.add_sequences("rep_seqs", fasta)
    with pytest.raises(MicrobiomeSuiteError, match="fasta asset expects"):
        dataset.add_sequences("bad_rep_seqs", not_fasta)


def test_dataset_rejects_mismatched_paired_fastq_formats(tmp_path: Path) -> None:
    r1 = tmp_path / "S1_R1.fastq.gz"
    r2 = tmp_path / "S1_R2.fq"
    r1.touch()
    r2.touch()

    dataset = MicrobiomeDataset()

    with pytest.raises(MicrobiomeSuiteError, match="same format"):
        dataset.add_reads("S1", r1=r1, r2=r2)


def test_dataset_public_api_exports() -> None:
    assert api.MicrobiomeDataset is MicrobiomeDataset
    assert api.DataAsset is DataAsset
    assert "MicrobiomeDataset" in api.__all__
    assert "DataAsset" in api.__all__
