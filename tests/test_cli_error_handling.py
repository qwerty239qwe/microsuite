from __future__ import annotations

import importlib

import pytest

from microsuite._errors import MicrobiomeSuiteError

cli_app = importlib.import_module("microsuite.cli.app")


def _raise_known_error() -> None:
    raise MicrobiomeSuiteError("expected failure")


def test_main_prints_known_error_without_traceback(monkeypatch, capsys) -> None:
    monkeypatch.setattr(cli_app, "app", _raise_known_error)
    monkeypatch.setattr(cli_app, "_DEBUG_ENABLED", False)

    with pytest.raises(SystemExit) as caught:
        cli_app.main()

    assert caught.value.code == 1
    stderr = capsys.readouterr().err
    assert "Error: expected failure" in stderr
    assert "Traceback" not in stderr
    assert "click.exceptions.Exit" not in stderr


def test_main_debug_reraises_known_error(monkeypatch) -> None:
    monkeypatch.setattr(cli_app, "app", _raise_known_error)
    monkeypatch.setattr(cli_app, "_DEBUG_ENABLED", True)

    with pytest.raises(MicrobiomeSuiteError, match="expected failure"):
        cli_app.main()
