from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata import write_resolved_config
from microsuite.metadata.schemas import RESOLVED_CONFIG_VERSION, published_schema_path
from microsuite.metadata.validate import validate


def test_write_resolved_config_redacts_and_validates(tmp_path: Path) -> None:
    cfg = {
        "project": {"accession": "ERP120510"},
        "dada2": {"max_ee_f": 2, "trunc_len_f": 240},
        "auth_token": "supersecret",
        "nested": {"api_key": "abc123"},
    }
    path = write_resolved_config(tmp_path, cfg)
    assert path.name == "resolved_config.json"
    data = json.loads(path.read_text())
    assert data["schema_version"] == RESOLVED_CONFIG_VERSION
    assert data["producer"]["name"] == "microsuite"
    assert data["config"]["dada2"]["max_ee_f"] == 2  # defaults preserved
    assert data["config"]["auth_token"] == "***"
    assert data["config"]["nested"]["api_key"] == "***"
    assert validate(data, RESOLVED_CONFIG_VERSION) == []
    assert not list(tmp_path.glob("*.tmp*"))


def test_write_resolved_config_custom_name(tmp_path: Path) -> None:
    path = write_resolved_config(tmp_path, {"a": 1}, name="snapshot.json")
    assert path.name == "snapshot.json"
    assert json.loads(path.read_text())["config"] == {"a": 1}


def test_write_resolved_config_rejects_non_mapping_config(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="invalid resolved-config"):
        write_resolved_config(tmp_path, [1, 2, 3])  # type: ignore[arg-type]


def test_resolved_config_published_json_schema_parity(tmp_path: Path) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    path = write_resolved_config(tmp_path, {"project": {"accession": "X"}})
    payload = json.loads(path.read_text())
    schema = json.loads(published_schema_path(RESOLVED_CONFIG_VERSION).read_text())
    cls = jsonschema.validators.validator_for(schema)
    cls.check_schema(schema)
    cls(schema, format_checker=cls.FORMAT_CHECKER).validate(payload)
    # a broken envelope is rejected by both validators
    payload["schema_version"] = "resolved-config.v2"
    assert validate(payload, RESOLVED_CONFIG_VERSION) != []
    assert not cls(schema).is_valid(payload)
