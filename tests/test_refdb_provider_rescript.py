from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.providers import get_provider
from microsuite.refdb.providers import rescript as _rescript  # noqa: F401  (force registration)
from microsuite.refdb.spec import RefDbSpec


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
        provider.fetch(RefDbSpec(name="silva", version="138.1", provider="rescript"), out_dir=tmp_path)
