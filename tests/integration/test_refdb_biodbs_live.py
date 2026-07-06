from __future__ import annotations

import os
from pathlib import Path

import pytest

pytest.importorskip("biodbs")

from microsuite.refdb.service import fetch_refdb  # noqa: E402
from microsuite.refdb.spec import RefDbSpec  # noqa: E402

pytestmark = [
    pytest.mark.skipif(
        os.environ.get("MICROSUITE_RUN_EXTERNAL_INTEGRATION") != "1",
        reason="set MICROSUITE_RUN_EXTERNAL_INTEGRATION=1 to run external-tool integration tests",
    ),
]


def _count_fasta_records(path: Path) -> int:
    return sum(1 for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">"))


def _fetch_and_check(monkeypatch, tmp_path: Path, name: str, version: str) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    artifact = fetch_refdb(RefDbSpec(name=name, version=version), "vsearch")

    assert artifact.path.exists()
    assert _count_fasta_records(artifact.path) >= 1


# --- light DBs (~11 MB each): full end-to-end, safe to run routinely ---


def test_homd_live_vsearch(monkeypatch, tmp_path: Path) -> None:
    _fetch_and_check(monkeypatch, tmp_path, "homd", "15.22")


def test_greengenes_live_vsearch(monkeypatch, tmp_path: Path) -> None:
    _fetch_and_check(monkeypatch, tmp_path, "greengenes", "2022.7-rc1")


# --- heavy DBs (39-200 MB each): gated the same way, but not run routinely ---


def test_pr2_live_vsearch(monkeypatch, tmp_path: Path) -> None:
    _fetch_and_check(monkeypatch, tmp_path, "pr2", "5.1.1")


def test_silva_live_vsearch(monkeypatch, tmp_path: Path) -> None:
    _fetch_and_check(monkeypatch, tmp_path, "silva", "138.2")


def test_gtdb_live_vsearch(monkeypatch, tmp_path: Path) -> None:
    _fetch_and_check(monkeypatch, tmp_path, "gtdb", "latest")


def test_unite_live_vsearch(monkeypatch, tmp_path: Path) -> None:
    _fetch_and_check(monkeypatch, tmp_path, "unite", "2020-02-20")
