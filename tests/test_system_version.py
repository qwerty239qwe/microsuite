from __future__ import annotations

import json

from typer.testing import CliRunner

from microsuite.cli.app import app
from microsuite.system import version as version_module


def test_version_info_has_stable_fields(monkeypatch) -> None:
    monkeypatch.setattr(version_module, "package_version", lambda: "0.2.0.dev0")
    monkeypatch.setattr(version_module, "source_commit", lambda: ("editable", "abc123"))

    payload = version_module.version_info()

    assert payload["name"] == "microsuite"
    assert payload["version"] == "0.2.0.dev0"
    assert payload["source"] == "editable"
    assert payload["commit"] == "abc123"
    assert payload["python"]


def test_version_cli_json() -> None:
    result = CliRunner().invoke(app, ["version", "--json"])

    assert result.exit_code == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["name"] == "microsuite"
    assert payload["version"]
