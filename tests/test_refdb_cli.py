# tests/test_refdb_cli.py
from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from microsuite.cli.app import app

FIXTURE = Path("src/microsuite/data/fixtures/refdb_mock")
runner = CliRunner()


def test_refdb_fetch_prints_artifact_path(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.setenv("MICROSUITE_REFDB_DIR", str(tmp_path / "cache"))

    class _FakeBiodbs:
        _SRC = {
            "HOMD_16S_rRNA_RefSeq_V16.03.fasta": FIXTURE / "source_a.fasta",
            "HOMD_16S_rRNA_RefSeq_V16.03.qiime.taxonomy": FIXTURE / "source_a.tax.tsv",
        }

        def homd_download_file(self, path_or_url, dest, overwrite=False):
            name = Path(path_or_url).name
            target = Path(dest) / name
            target.write_text(self._SRC[name].read_text(encoding="utf-8"), encoding="utf-8")
            return target

    monkeypatch.setattr(
        "microsuite.refdb.providers.biodbs._load_biodbs", lambda: _FakeBiodbs()
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
