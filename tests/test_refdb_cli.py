# tests/test_refdb_cli.py
from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from microsuite.cli.app import app

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")
runner = CliRunner()


def test_refdb_fetch_prints_artifact_path(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))

    def fake_fetch(name, version, out_dir):
        return (str(FIXTURE / "source_a.fasta"), str(FIXTURE / "source_a.tax.tsv"))

    monkeypatch.setattr(
        "microsuite.refdb.providers.biodbs._load_biodbs_fetch", lambda: fake_fetch
    )
    result = runner.invoke(
        app, ["refdb", "fetch", "homd", "--version", "15.22", "--build", "vsearch"]
    )
    assert result.exit_code == 0, result.output
    assert "reference.fasta" in result.output


def test_refdb_fetch_rejects_bad_build(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))
    result = runner.invoke(
        app, ["refdb", "fetch", "homd", "--version", "15.22", "--build", "bowtie"]
    )
    assert result.exit_code != 0
