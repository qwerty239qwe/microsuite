from __future__ import annotations

import json

from typer.testing import CliRunner

from microsuite.cli.app import app
from microsuite.system.capabilities import capability_payload, require_capabilities


def test_capability_payload_contract() -> None:
    payload = capability_payload()

    assert payload["schema_version"] == "microsuite-capabilities.v1"
    assert payload["producer"]["name"] == "microsuite"
    capability = payload["capabilities"]["diversity.beta_significance.permdisp.native"]
    assert capability == {"available": True, "api": 1}


def test_require_capabilities_reports_missing() -> None:
    available, missing = require_capabilities(
        ["diversity.beta_significance.permdisp.native", "not.available"]
    )

    assert available == ["diversity.beta_significance.permdisp.native"]
    assert missing == ["not.available"]


def test_capabilities_cli_json() -> None:
    result = CliRunner().invoke(app, ["capabilities", "--json"])

    assert result.exit_code == 0, result.stdout
    payload = json.loads(result.stdout)
    assert payload["schema_version"] == "microsuite-capabilities.v1"
