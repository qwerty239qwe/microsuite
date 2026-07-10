from __future__ import annotations

import gzip
import io
import tarfile
import zipfile
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import biodbs as _biodbs  # noqa: F401  (registration)
from microsuite.refdb.providers import get_provider
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
    ids = [
        line[1:].strip() for line in raw.sequences.read_text().splitlines() if line.startswith(">")
    ]
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
        body = (
            b">AY855839.1.1390 Bacteria;Firmicutes;Bacilli\nACGT\n"
            b">FJ12.1.1500 Bacteria;Bacteroidetes\nTTTT\n"
        )
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
    ids = [
        line[1:].strip() for line in raw.sequences.read_text().splitlines() if line.startswith(">")
    ]
    assert ids == ["AY855839.1.1390", "FJ12.1.1500"]  # id only, lineage stripped from header
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    assert tax["AY855839.1.1390"] == "Bacteria;Firmicutes;Bacilli"


def test_gtdb_adapter_joins_ssu_to_taxonomy_tilde_edge_case(tmp_path: Path, monkeypatch) -> None:
    # Edge case: some SSU headers carry a '~contig' suffix after the genome accession.
    # Not observed in live GTDB data (see test below for the real shape), but the adapter's
    # split("~", 1)[0] must still resolve the lineage correctly when it does occur.
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeSilvaGtdb(tmp_path))
    raw = get_provider("biodbs").fetch(RefDbSpec(name="gtdb", version="latest"), out_dir=tmp_path)
    ids = [
        line[1:].split()[0]
        for line in raw.sequences.read_text().splitlines()
        if line.startswith(">")
    ]
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
    ids = [
        line[1:].split()[0]
        for line in raw.sequences.read_text().splitlines()
        if line.startswith(">")
    ]
    assert ids == ["RS_GCF_1", "RS_GCF_2"]
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    assert tax["RS_GCF_1"] == "d__Bacteria;p__Firmicutes"
    assert tax["RS_GCF_2"] == "d__Bacteria;p__Actinobacteria"


def _write_qza(zip_path: Path, member: str, content: str) -> None:
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr(f"5cc3556c-5f05-4cac-8d3b-8727a0d3b705/{member}", content)


class _FakeGreenGenes:
    """Mimics greengenes_download_file: writes tiny real .qza (zip) fixtures,
    matching the live-verified QIIME2 artifact layout (<uuid>/data/<name>)."""

    def greengenes_download_file(self, path, dest, overwrite=False) -> Path:
        name = Path(path).name
        target = Path(dest) / name
        if name.endswith(".fna.qza"):
            _write_qza(target, "data/dna-sequences.fasta", ">ggA\nACGT\n>ggB\nTGCA\n")
        elif name.endswith(".tax.qza"):
            _write_qza(
                target,
                "data/taxonomy.tsv",
                "Feature ID\tTaxon\n"
                "ggA\td__Bacteria; p__Firmicutes\n"
                "ggB\td__Bacteria; p__Bacteroidota\n",
            )
        else:
            raise AssertionError(f"unexpected greengenes path: {path}")
        return target


def test_greengenes_adapter_extracts_qza_members(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeGreenGenes())
    raw = get_provider("biodbs").fetch(
        RefDbSpec(name="greengenes", version="2022.7-rc1"), out_dir=tmp_path
    )
    ids = [
        line[1:].strip() for line in raw.sequences.read_text().splitlines() if line.startswith(">")
    ]
    assert ids == ["ggA", "ggB"]
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    assert tax["ggA"] == "d__Bacteria; p__Firmicutes"
    assert tax["ggB"] == "d__Bacteria; p__Bacteroidota"
    # header row must be stripped, not treated as a record
    assert "Feature ID" not in tax


def test_greengenes_rejects_unsupported_version(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeGreenGenes())
    with pytest.raises(MicrobiomeSuiteError):
        get_provider("biodbs").fetch(
            RefDbSpec(name="greengenes", version="2024.09"), out_dir=tmp_path
        )


class _FakeUnite:
    """Mimics unite_download: writes a tiny real .tgz fixture that mirrors the
    live-verified UNITE archive layout — multiple fasta/taxonomy candidates
    (dynamic + 97/99 thresholds, plus a developer/ duplicate set) so the
    adapter must disambiguate rather than assume a single match."""

    def unite_download(self, version, dest, taxon_group="fungi", singletons=False) -> Path:
        target = Path(dest) / "unite.tgz"
        root = "sh_qiime_release_04.02.2020"
        with tarfile.open(target, "w:gz") as tf:

            def add(name: str, content: bytes) -> None:
                info = tarfile.TarInfo(name=name)
                info.size = len(content)
                tf.addfile(info, io.BytesIO(content))

            add(
                f"{root}/sh_refs_qiime_ver8_dynamic_04.02.2020.fasta",
                b">SH001.08FU_UN001_refs\nACGT\n>SH002.08FU_UN002_refs\nTGCA\n",
            )
            add(
                f"{root}/sh_taxonomy_qiime_ver8_dynamic_04.02.2020.txt",
                b"SH001.08FU_UN001_refs\tk__Fungi;p__Ascomycota\n"
                b"SH002.08FU_UN002_refs\tk__Fungi;p__Basidiomycota\n",
            )
            # decoy 97%-threshold pair (must NOT be picked over the dynamic one)
            add(
                f"{root}/sh_refs_qiime_ver8_97_04.02.2020.fasta",
                b">decoyA\nAAAA\n",
            )
            add(
                f"{root}/sh_taxonomy_qiime_ver8_97_04.02.2020.txt",
                b"decoyA\tk__Fungi;p__Decoyphyta\n",
            )
            # developer/ duplicate of the dynamic pair (must NOT be picked either)
            add(
                f"{root}/developer/sh_refs_qiime_ver8_dynamic_04.02.2020_dev.fasta",
                b">devA\nCCCC\n",
            )
            add(
                f"{root}/developer/sh_taxonomy_qiime_ver8_dynamic_04.02.2020_dev.txt",
                b"devA\tk__Fungi;p__Devphyta\n",
            )
        return target


def test_unite_adapter_picks_dynamic_pair_from_tgz(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakeUnite())
    raw = get_provider("biodbs").fetch(
        RefDbSpec(name="unite", version="2020-02-20"), out_dir=tmp_path
    )
    ids = [
        line[1:].strip() for line in raw.sequences.read_text().splitlines() if line.startswith(">")
    ]
    assert ids == ["SH001.08FU_UN001_refs", "SH002.08FU_UN002_refs"]
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    assert tax["SH001.08FU_UN001_refs"] == "k__Fungi;p__Ascomycota"
    assert tax["SH002.08FU_UN002_refs"] == "k__Fungi;p__Basidiomycota"
    assert "decoyA" not in tax
    assert "devA" not in tax


class _FakePr2:
    """Mimics pr2_download_asset: writes tiny real .gz mothur fixtures, matching
    the live-verified id-first `accession<TAB>lineage;` (trailing ';') shape."""

    def pr2_download_asset(self, name, dest, tag=None) -> Path:
        target = Path(dest) / name
        if name.endswith("mothur.fasta.gz"):
            body = b">pr2A.1.1_U\nACGT\n>pr2B.1.1_U\nTGCA\n"
        elif name.endswith("mothur.tax.gz"):
            body = (
                b"pr2A.1.1_U\tEukaryota;TSAR;Alveolata;\n"
                b"pr2B.1.1_U\tEukaryota;Obazoa;Opisthokonta;\n"
            )
        else:
            raise AssertionError(f"unexpected pr2 asset: {name}")
        target.write_bytes(gzip.compress(body))
        return target


def test_pr2_adapter_parses_mothur_files(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr(_biodbs, "_load_biodbs", lambda: _FakePr2())
    raw = get_provider("biodbs").fetch(RefDbSpec(name="pr2", version="5.1.1"), out_dir=tmp_path)
    ids = [
        line[1:].strip() for line in raw.sequences.read_text().splitlines() if line.startswith(">")
    ]
    assert ids == ["pr2A.1.1_U", "pr2B.1.1_U"]
    tax = dict(r.split("\t", 1) for r in raw.taxonomy.read_text().splitlines() if r.strip())
    assert tax["pr2A.1.1_U"] == "Eukaryota;TSAR;Alveolata"
    assert tax["pr2B.1.1_U"] == "Eukaryota;Obazoa;Opisthokonta"
