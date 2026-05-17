from __future__ import annotations

import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.diffab import ancombc
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_tsv
from microsuite.methods.diff_abundance import diff_abundance

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_table(tmp_path: Path) -> Path:
    table = tmp_path / "table.h5ad"
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    write_h5ad(adata, table)
    return table


def test_ancombc_r_script_is_external_asset() -> None:
    script = ROOT / "scripts" / "r" / "ancombc.R"

    assert script.exists()
    text = script.read_text(encoding="utf-8")
    assert "ANCOMBC" in text
    assert "commandArgs" in text


def test_diff_abundance_ancombc_missing_rscript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="Rscript"):
        diff_abundance(
            backend="ancombc",
            table=table,
            group="treatment",
            output=tmp_path / "diff.tsv",
        )


def test_ancombc_invokes_external_script_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    commands: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript")

    def fake_run(command: list[str], **kwargs: object) -> object:
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(ancombc.subprocess, "run", fake_run)

    diff_abundance(
        backend="ancombc",
        table=table,
        group="treatment",
        output=tmp_path / "diff.tsv",
    )

    assert commands
    command = commands[0]
    assert command[0] == "Rscript"
    assert command[1] == str(ROOT / "scripts" / "r" / "ancombc.R")
    assert command[-2:] == ["treatment", str(tmp_path / "diff.tsv")]


def test_diff_abundance_planned_backend_message(tmp_path: Path) -> None:
    table = fixture_table(tmp_path)

    with pytest.raises(MicrobiomeSuiteError, match="not implemented"):
        diff_abundance(
            backend="aldex2",
            table=table,
            group="treatment",
            output=tmp_path / "diff.tsv",
        )


def test_cli_diff_abundance_help_and_missing_rscript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: None)
    runner = CliRunner()

    help_result = runner.invoke(app, ["diff_abundance", "--help"])
    assert help_result.exit_code == 0
    assert "--backend" in help_result.stdout

    result = runner.invoke(
        app,
        [
            "diff_abundance",
            "--backend",
            "ancombc",
            "--table",
            str(table),
            "--group",
            "treatment",
            "-o",
            str(tmp_path / "diff.tsv"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
    assert "Rscript" in str(result.exception)


def test_legacy_diffab_ancombc_command_still_reports_rscript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: None)

    result = CliRunner().invoke(
        app,
        [
            "diffab",
            "ancombc",
            str(table),
            "--group",
            "treatment",
            "-o",
            str(tmp_path / "diff.tsv"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
    assert "Rscript" in str(result.exception)
