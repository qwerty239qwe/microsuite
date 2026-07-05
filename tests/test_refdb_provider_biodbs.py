from __future__ import annotations

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
