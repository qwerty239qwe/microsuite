from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider
from microsuite.refdb.providers import biodbs as _biodbs  # noqa: F401  (force registration)
from microsuite.refdb.spec import RefDbSpec

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")


def test_biodbs_is_default_provider() -> None:
    assert RefDbSpec(name="homd", version="15.22").provider == "biodbs"


def test_biodbs_fetch_uses_upstream_api(tmp_path: Path, monkeypatch) -> None:
    def fake_fetch(name, version, out_dir):
        return (str(FIXTURE / "source_a.fasta"), str(FIXTURE / "source_a.tax.tsv"))

    monkeypatch.setattr(_biodbs, "_load_biodbs_fetch", lambda: fake_fetch)
    provider = get_provider("biodbs")
    raw = provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)
    assert raw.sequences.name == "source_a.fasta"
    assert raw.taxonomy.name == "source_a.tax.tsv"


def test_biodbs_missing_dependency_raises(tmp_path: Path, monkeypatch) -> None:
    def boom():
        raise ImportError("no biodbs")

    monkeypatch.setattr(_biodbs, "_load_biodbs_fetch", boom)
    provider = get_provider("biodbs")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(RefDbSpec(name="homd", version="15.22"), out_dir=tmp_path)
