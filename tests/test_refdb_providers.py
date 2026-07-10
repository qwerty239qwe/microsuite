from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider, register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


class FakeProvider(RefDbProvider):
    name = "fake"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        return RawRefDb(
            sequences=FIXTURE / "source_a.fasta",
            taxonomy=FIXTURE / "source_a.tax.tsv",
        )


def test_register_and_get_provider() -> None:
    register_provider(FakeProvider())
    assert isinstance(get_provider("fake"), FakeProvider)


def test_get_unknown_provider_raises() -> None:
    with pytest.raises(MicrobiomeSuiteError):
        get_provider("does-not-exist")


def test_default_build_delegates_to_build_artifact(tmp_path: Path) -> None:
    provider = FakeProvider()
    raw = provider.fetch(RefDbSpec(name="x", version="1"), out_dir=tmp_path)
    art = provider.build(raw, "vsearch", out_dir=tmp_path)
    assert art.build_target == "vsearch"
    assert art.path.exists()


class QzaProvider(RefDbProvider):
    """Provider whose fetch() already yields a packaged QIIME2 .qza."""

    name = "qza-fake"

    def __init__(self, qza_path: Path) -> None:
        self._qza_path = qza_path

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        return RawRefDb(
            sequences=self._qza_path,
            taxonomy=FIXTURE / "source_a.tax.tsv",
            qza=self._qza_path,
        )


def _make_qza(tmp_path: Path) -> Path:
    qza_path = tmp_path / "packaged.qza"
    qza_path.write_bytes(b"fake-zip-bytes")
    return qza_path


def test_build_qiime2_short_circuits_on_existing_qza(tmp_path: Path, monkeypatch) -> None:
    qza_path = _make_qza(tmp_path)
    provider = QzaProvider(qza_path)
    raw = provider.fetch(RefDbSpec(name="x", version="1"), out_dir=tmp_path)

    called = False

    def fake_build_artifact(*args, **kwargs):
        nonlocal called
        called = True
        raise AssertionError("build_artifact must not be called for qza short-circuit")

    monkeypatch.setattr("microsuite.refdb.providers._base.build_artifact", fake_build_artifact)

    art = provider.build(raw, "qiime2", out_dir=tmp_path)

    assert not called
    assert art.build_target == "qiime2"
    assert art.path == qza_path
    assert len(art.checksum) == 64
    int(art.checksum, 16)  # valid hex


def test_build_vsearch_raises_on_existing_qza(tmp_path: Path) -> None:
    qza_path = _make_qza(tmp_path)
    provider = QzaProvider(qza_path)
    raw = provider.fetch(RefDbSpec(name="x", version="1"), out_dir=tmp_path)

    with pytest.raises(MicrobiomeSuiteError):
        provider.build(raw, "vsearch", out_dir=tmp_path)


def test_build_blast_raises_on_existing_qza(tmp_path: Path) -> None:
    qza_path = _make_qza(tmp_path)
    provider = QzaProvider(qza_path)
    raw = provider.fetch(RefDbSpec(name="x", version="1"), out_dir=tmp_path)

    with pytest.raises(MicrobiomeSuiteError):
        provider.build(raw, "blast", out_dir=tmp_path)


def test_build_vsearch_without_qza_still_delegates(tmp_path: Path) -> None:
    provider = FakeProvider()
    raw = provider.fetch(RefDbSpec(name="x", version="1"), out_dir=tmp_path)
    assert raw.qza is None
    art = provider.build(raw, "vsearch", out_dir=tmp_path)
    assert art.build_target == "vsearch"
    assert art.path.exists()
