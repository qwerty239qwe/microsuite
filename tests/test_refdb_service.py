# tests/test_refdb_service.py
from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import register_provider
from microsuite.refdb.providers._base import RefDbProvider
from microsuite.refdb.registry import RefDbRegistry
from microsuite.refdb.service import fetch_refdb, resolve_classifier
from microsuite.refdb.spec import RawRefDb, RefDbSource, RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


class CountingProvider(RefDbProvider):
    name = "counting"

    def __init__(self) -> None:
        self.fetch_calls = 0

    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb:
        self.fetch_calls += 1
        tag = spec.name if (FIXTURE / f"{spec.name}.fasta").exists() else "source_a"
        return RawRefDb(
            sequences=FIXTURE / f"{tag}.fasta",
            taxonomy=FIXTURE / f"{tag}.tax.tsv",
        )


def test_fetch_then_cache_skips_second_fetch(tmp_path: Path) -> None:
    provider = CountingProvider()
    register_provider(provider)
    reg = RefDbRegistry(tmp_path / "cache")
    spec = RefDbSpec(name="source_a", version="1", provider="counting", build_targets=("vsearch",))

    first = fetch_refdb(spec, "vsearch", registry=reg)
    second = fetch_refdb(spec, "vsearch", registry=reg)

    assert first.path == second.path
    assert provider.fetch_calls == 1  # second call served from cache


def test_fetch_merges_multiple_sources(tmp_path: Path) -> None:
    register_provider(CountingProvider())
    reg = RefDbRegistry(tmp_path / "cache")
    spec = RefDbSpec(
        name="fomc-combined",
        version="20221029",
        provider="counting",
        sources=(RefDbSource("source_a", "1"), RefDbSource("source_b", "1")),
    )
    art = fetch_refdb(spec, "vsearch", registry=reg)
    ids = [line[1:].strip() for line in art.path.read_text().splitlines() if line.startswith(">")]
    assert ids == ["seq1", "seq2", "seq3"]


def test_resolve_classifier_raw_path_passthrough(tmp_path: Path) -> None:
    raw = tmp_path / "my.qza"
    raw.write_text("x", encoding="utf-8")
    assert resolve_classifier(str(raw)) == raw


def test_resolve_classifier_registry_ref(tmp_path: Path) -> None:
    register_provider(CountingProvider())
    reg = RefDbRegistry(tmp_path / "cache")
    spec = RefDbSpec(name="source_a", version="1", provider="counting")
    art = fetch_refdb(spec, "vsearch", registry=reg)
    resolved = resolve_classifier("refdb:source_a@1", registry=reg)
    assert resolved == art.path


def test_resolve_classifier_unknown_ref_raises(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    with pytest.raises(MicrobiomeSuiteError):
        resolve_classifier("refdb:ghost@9", registry=reg)
