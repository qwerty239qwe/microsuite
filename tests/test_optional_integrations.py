from __future__ import annotations

import importlib
import shutil
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diffab.ancombc import run_ancombc
from microsuite.io.biom import read_biom
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def test_biom_missing_dependency_message() -> None:
    pytest.importorskip(
        "biom", reason="biom-format is installed; missing-dependency branch skipped"
    )
    pytest.skip("biom-format installed")


def test_biom_import_reports_missing_dependency_when_unavailable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_import_module = importlib.import_module

    def fake_import_module(name: str, package: str | None = None) -> object:
        if name == "biom":
            raise ImportError("blocked")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)
    with pytest.raises(MicrobiomeSuiteError, match="biom-format"):
        read_biom(FIXTURE / "table.biom", FIXTURE / "metadata.tsv")


def test_ancombc_missing_rscript_message(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    if shutil.which("Rscript") is not None:
        pytest.skip("Rscript is installed")
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    monkeypatch.setattr(shutil, "which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="Rscript"):
        run_ancombc(adata, group="treatment", output=tmp_path / "ancombc.tsv")
