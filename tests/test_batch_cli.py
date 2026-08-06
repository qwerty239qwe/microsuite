from __future__ import annotations

from pathlib import Path

from typer.testing import CliRunner

from microsuite.cli.app import app

runner = CliRunner()


def test_batch_correct_is_registered() -> None:
    result = runner.invoke(app, ["batch", "--help"])
    assert result.exit_code == 0
    assert "correct" in result.stdout


def test_options_reach_the_method(monkeypatch, tmp_path: Path) -> None:
    from microsuite.cli import batch_cmd

    captured: dict = {}
    monkeypatch.setattr(batch_cmd, "batch_correct", lambda **kw: captured.update(kw))
    table = tmp_path / "t.h5ad"
    table.write_bytes(b"")
    result = runner.invoke(
        app,
        [
            "batch",
            "correct",
            str(table),
            "--output",
            str(tmp_path / "out.h5ad"),
            "--batch-col",
            "run_id",
            "--covariates",
            "sex",
            "--covariates",
            "age",
            "--backend",
            "mmuphin",
            "--runtime",
            "docker",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert captured["batch"] == "run_id"
    assert captured["covariates"] == ["sex", "age"]
    assert captured["backend"] == "mmuphin"
    assert captured["runtime"] == "docker"


def test_default_backend_is_mmuphin(monkeypatch, tmp_path: Path) -> None:
    from microsuite.cli import batch_cmd

    captured: dict = {}
    monkeypatch.setattr(batch_cmd, "batch_correct", lambda **kw: captured.update(kw))
    table = tmp_path / "t.h5ad"
    table.write_bytes(b"")
    result = runner.invoke(
        app,
        [
            "batch",
            "correct",
            str(table),
            "--output",
            str(tmp_path / "o.h5ad"),
            "--batch-col",
            "run_id",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert captured["backend"] == "mmuphin"
