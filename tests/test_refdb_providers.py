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
