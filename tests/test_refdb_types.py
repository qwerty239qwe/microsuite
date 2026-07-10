from __future__ import annotations

from pathlib import Path

from microsuite.refdb.paths import VALID_BUILD_TARGETS, refdb_cache_dir
from microsuite.refdb.spec import RefDbSpec


def test_refdbspec_defaults() -> None:
    spec = RefDbSpec(name="homd", version="15.22")
    assert spec.provider == "biodbs"
    assert spec.target == "16S"
    assert spec.build_targets == ("vsearch",)
    assert spec.sources == ()


def test_refdbspec_is_frozen() -> None:
    spec = RefDbSpec(name="homd", version="15.22")
    try:
        spec.name = "other"  # type: ignore[misc]
    except AttributeError:
        return
    raise AssertionError("RefDbSpec should be frozen")


def test_valid_build_targets() -> None:
    assert VALID_BUILD_TARGETS == ("vsearch", "blast", "qiime2")


def test_cache_dir_honors_env(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    root = refdb_cache_dir()
    assert root == tmp_path / "cache"
    assert root.is_dir()


def test_cache_dir_default(monkeypatch) -> None:
    monkeypatch.delenv("MICROSUITE_REFDB_DIR", raising=False)
    root = refdb_cache_dir()
    assert root == Path.home() / ".cache" / "microsuite" / "refdb"
