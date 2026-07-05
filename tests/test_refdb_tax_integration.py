# tests/test_refdb_tax_integration.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.tax_classify import tax_classify
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.registry import RefDbRegistry
from microsuite.refdb.service import fetch_refdb
from microsuite.refdb.spec import RawRefDb, RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


class FixtureProvider(RefDbProvider):
    name = "fixture"

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        return RawRefDb(
            sequences=FIXTURE / "source_a.fasta",
            taxonomy=FIXTURE / "source_a.tax.tsv",
        )


def test_tax_classify_resolves_refdb_ref(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    register_provider(FixtureProvider())
    reg = RefDbRegistry(tmp_path / "cache")
    art = fetch_refdb(
        RefDbSpec(name="source_a", version="1", provider="fixture"), "vsearch", registry=reg
    )

    captured: dict[str, Path | None] = {}

    def fake_qiime(*, rep_seqs, classifier, output, threads, force, run_dir, timeout):
        captured["classifier"] = classifier

    monkeypatch.setattr("microsuite.methods.tax_classify.tax_classify_qiime2", fake_qiime)
    rep = tmp_path / "rep.qza"
    rep.write_text("x", encoding="utf-8")

    tax_classify(
        backend="qiime2",
        rep_seqs=rep,
        classifier="refdb:source_a@1",
        output=tmp_path / "out.qza",
        threads=1,
        force=True,
    )
    assert captured["classifier"] == art.path


def test_tax_classify_unknown_ref_raises(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    rep = tmp_path / "rep.qza"
    rep.write_text("x", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        tax_classify(
            backend="qiime2",
            rep_seqs=rep,
            classifier="refdb:ghost@9",
            output=tmp_path / "out.qza",
            threads=1,
            force=True,
        )
