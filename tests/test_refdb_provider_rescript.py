from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider
from microsuite.refdb.providers import rescript as _rescript  # noqa: F401  (force registration)
from microsuite.refdb.providers.rescript import RescriptProvider
from microsuite.refdb.spec import RawRefDb, RefDbSpec


def test_rescript_silva_builds_command(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: "qiime" if name == "qiime" else None)
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        for i, tok in enumerate(command):
            if tok == "--o-silva-sequences":
                Path(command[i + 1]).write_text("x", encoding="utf-8")
            if tok == "--o-silva-taxonomy":
                Path(command[i + 1]).write_text("x", encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)
    provider = get_provider("rescript")
    provider.fetch(RefDbSpec(name="silva", version="138.1", provider="rescript"), out_dir=tmp_path)

    assert calls[0][:3] == ["qiime", "rescript", "get-silva-data"]


def test_rescript_requires_qiime(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setattr("shutil.which", lambda name: None)
    provider = get_provider("rescript")
    with pytest.raises(MicrobiomeSuiteError):
        provider.fetch(
            RefDbSpec(name="silva", version="138.1", provider="rescript"), out_dir=tmp_path
        )


def test_rescript_build_qiime2_short_circuits_on_qza(tmp_path: Path) -> None:
    seqs = tmp_path / "x.qza"
    tax = tmp_path / "t.qza"
    seqs.write_bytes(b"fake-seqs-zip")
    tax.write_bytes(b"fake-tax-zip")
    raw = RawRefDb(sequences=seqs, taxonomy=tax, qza=seqs)

    art = RescriptProvider().build(raw, "qiime2", tmp_path)

    assert art.path == seqs
    assert art.build_target == "qiime2"


def test_rescript_build_vsearch_raises_on_qza(tmp_path: Path) -> None:
    seqs = tmp_path / "x.qza"
    tax = tmp_path / "t.qza"
    seqs.write_bytes(b"fake-seqs-zip")
    tax.write_bytes(b"fake-tax-zip")
    raw = RawRefDb(sequences=seqs, taxonomy=tax, qza=seqs)

    with pytest.raises(MicrobiomeSuiteError):
        RescriptProvider().build(raw, "vsearch", tmp_path)
