from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider
from microsuite.refdb.providers import biodbs as _biodbs  # noqa: F401  (registration)
from microsuite.refdb.spec import RefDbSpec


class _FakeBiodbs:
    """Mimics biodbs.homd_download_file: dest is a dir, basename taken from the URL path.
    Writes fixture content for the two HOMD files the adapter requests."""

    _CONTENT = {
        "HOMD_16S_rRNA_RefSeq_V16.03.fasta": ">seqA\nACGT\n>seqB\nTGCA\n",
        "HOMD_16S_rRNA_RefSeq_V16.03.qiime.taxonomy": "seqA\tk__B;s__x\nseqB\tk__B;s__y\n",
    }

    def homd_download_file(self, path_or_url, dest, overwrite=False) -> Path:
        name = Path(path_or_url).name
        target = Path(dest) / name
        target.write_text(self._CONTENT[name], encoding="utf-8")
        return target


def test_biodbs_is_default_provider() -> None:
    assert RefDbSpec(name="homd", version="15.22").provider == "biodbs"


def test_homd_adapter_produces_seqs_and_taxonomy(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeBiodbs())
    provider = get_provider("biodbs")
    raw = provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)
    assert raw.sequences.exists()
    ids = [l[1:].strip() for l in raw.sequences.read_text().splitlines() if l.startswith(">")]
    assert ids == ["seqA", "seqB"]
    tax_first_col = [r.split("\t")[0] for r in raw.taxonomy.read_text().splitlines() if r.strip()]
    assert tax_first_col == ["seqA", "seqB"]


def test_unknown_db_raises(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeBiodbs())
    provider = get_provider("biodbs")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(RefDbSpec(name="not-a-db", version="1"), out_dir=tmp_path)


def test_missing_biodbs_raises(tmp_path: Path, monkeypatch) -> None:
    def boom():
        raise MicrobiomeSuiteError("biodbs not installed")

    monkeypatch.setattr(_biodbs, "_load_biodbs", boom)
    provider = get_provider("biodbs")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)


class _FakeSilvaGtdb:
    def __init__(self, tmp: Path) -> None:
        self._tmp = tmp

    def silva_download_file(self, url, dest, overwrite=False) -> Path:
        target = Path(dest) / Path(url).name  # ...SSURef_NR99_tax_silva.fasta.gz
        body = b">AY855839.1.1390 Bacteria;Firmicutes;Bacilli\nACGT\n>FJ12.1.1500 Bacteria;Bacteroidetes\nTTTT\n"
        target.write_bytes(gzip.compress(body))
        return target

    def gtdb_download_file(self, path, dest, overwrite=False) -> Path:
        target = Path(dest) / Path(path).name  # bac120_ssu_reps.fna.gz
        body = b">RS_GCF_1~ctg1 desc\nACGT\n>RS_GCF_2~ctg9 desc\nGGGG\n"
        target.write_bytes(gzip.compress(body))
        return target

    def gtdb_download_taxonomy(
        self, domain, dest, release="latest", compressed=True, overwrite=False
    ) -> Path:
        target = Path(dest) / f"{domain}_taxonomy.tsv.gz"
        body = b"RS_GCF_1\td__Bacteria;p__Firmicutes\nRS_GCF_2\td__Bacteria;p__Actinobacteria\n"
        target.write_bytes(gzip.compress(body))
        return target


def test_silva_adapter_parses_headers(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeSilvaGtdb(tmp_path))
    raw = get_provider("biodbs").fetch(RefDbSpec(name="silva", version="138.2"), out_dir=tmp_path)
    ids = [l[1:].strip() for l in raw.sequences.read_text().splitlines() if l.startswith(">")]
    assert ids == ["AY855839.1.1390", "FJ12.1.1500"]  # id only, lineage stripped from header
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    assert tax["AY855839.1.1390"] == "Bacteria;Firmicutes;Bacilli"


def test_gtdb_adapter_joins_ssu_to_taxonomy_tilde_edge_case(tmp_path: Path, monkeypatch) -> None:
    # Edge case: some SSU headers carry a '~contig' suffix after the genome accession.
    # Not observed in live GTDB data (see test below for the real shape), but the adapter's
    # split("~", 1)[0] must still resolve the lineage correctly when it does occur.
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeSilvaGtdb(tmp_path))
    raw = get_provider("biodbs").fetch(RefDbSpec(name="gtdb", version="latest"), out_dir=tmp_path)
    ids = [l[1:].split()[0] for l in raw.sequences.read_text().splitlines() if l.startswith(">")]
    assert ids == ["RS_GCF_1~ctg1", "RS_GCF_2~ctg9"]
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    # taxonomy keyed by the SSU record id, lineage looked up via genome accession before '~'
    assert tax["RS_GCF_1~ctg1"] == "d__Bacteria;p__Firmicutes"
    assert tax["RS_GCF_2~ctg9"] == "d__Bacteria;p__Actinobacteria"


class _FakeSilvaGtdbBareAccession:
    """Mimics REAL GTDB SSU headers: the record id is the bare genome accession with no
    '~contig' suffix, matching the taxonomy TSV directly. Verified against live GTDB output
    in docs/superpowers/specs/biodbs-v040-probe.md (e.g. '>RS_GCF_031457235.1 d__Bacteria;...')."""

    def __init__(self, tmp: Path) -> None:
        self._tmp = tmp

    def gtdb_download_file(self, path, dest, overwrite=False) -> Path:
        target = Path(dest) / Path(path).name
        body = b">RS_GCF_1 desc\nACGT\n>RS_GCF_2 desc\nGGGG\n"
        target.write_bytes(gzip.compress(body))
        return target

    def gtdb_download_taxonomy(
        self, domain, dest, release="latest", compressed=True, overwrite=False
    ) -> Path:
        target = Path(dest) / f"{domain}_taxonomy.tsv.gz"
        body = b"RS_GCF_1\td__Bacteria;p__Firmicutes\nRS_GCF_2\td__Bacteria;p__Actinobacteria\n"
        target.write_bytes(gzip.compress(body))
        return target


def test_gtdb_adapter_joins_ssu_to_taxonomy_real_bare_accession(
    tmp_path: Path, monkeypatch
) -> None:
    # Production-shaped case: bare-accession SSU headers with no '~' separator, as confirmed
    # by live GTDB verification. This is the format that actually ships.
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeSilvaGtdbBareAccession(tmp_path))
    raw = get_provider("biodbs").fetch(RefDbSpec(name="gtdb", version="latest"), out_dir=tmp_path)
    ids = [l[1:].split()[0] for l in raw.sequences.read_text().splitlines() if l.startswith(">")]
    assert ids == ["RS_GCF_1", "RS_GCF_2"]
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    assert tax["RS_GCF_1"] == "d__Bacteria;p__Firmicutes"
    assert tax["RS_GCF_2"] == "d__Bacteria;p__Actinobacteria"
