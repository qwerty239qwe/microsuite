from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.registry import RefDbRegistry, sha256_file
from microsuite.refdb.spec import BuiltArtifact


def _artifact(tmp_path: Path, text: str = "ACGT") -> BuiltArtifact:
    path = tmp_path / "db.fasta"
    path.write_text(text, encoding="utf-8")
    return BuiltArtifact(path=path, build_target="vsearch", checksum=sha256_file(path))


def test_record_then_resolve_round_trip(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    art = _artifact(tmp_path)
    reg.record("homd", "15.22", art, provider="biodbs")

    resolved = reg.resolve("homd", "15.22", "vsearch")
    assert resolved is not None
    assert resolved.path == art.path
    assert resolved.checksum == art.checksum


def test_resolve_missing_returns_none(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    assert reg.resolve("nope", "1", "vsearch") is None


def test_resolve_checksum_mismatch_returns_none(tmp_path: Path) -> None:
    reg = RefDbRegistry(tmp_path / "cache")
    art = _artifact(tmp_path)
    reg.record("homd", "15.22", art, provider="biodbs")
    art.path.write_text("MUTATED", encoding="utf-8")  # invalidate
    assert reg.resolve("homd", "15.22", "vsearch") is None


def test_corrupt_manifest_raises(tmp_path: Path) -> None:
    root = tmp_path / "cache"
    root.mkdir(parents=True)
    (root / "manifest.json").write_text("{ not json", encoding="utf-8")
    reg = RefDbRegistry(root)
    with pytest.raises(MicrobiomeSuiteError):
        reg.resolve("homd", "15.22", "vsearch")
