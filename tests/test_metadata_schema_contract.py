from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite.metadata.schemas import published_schema_path
from microsuite.metadata.validate import validate, validate_stage_result

FIXTURES = Path(__file__).parent / "fixtures" / "stage_result"
VALID = sorted((FIXTURES / "valid").glob("*.json"))
INVALID = sorted((FIXTURES / "invalid").glob("*.json"))


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_published_schema_is_valid_json() -> None:
    schema = json.loads(published_schema_path().read_text(encoding="utf-8"))
    assert schema["$id"].endswith("stage-result.v1.schema.json")
    assert schema["properties"]["schema_version"]["const"] == "stage-result.v1"


@pytest.mark.parametrize("path", VALID, ids=[p.stem for p in VALID])
def test_valid_fixtures_pass_python_validator(path: Path) -> None:
    assert validate_stage_result(_load(path)) == []


@pytest.mark.parametrize("path", INVALID, ids=[p.stem for p in INVALID])
def test_invalid_fixtures_fail_python_validator(path: Path) -> None:
    assert validate(_load(path)) != []


def _jsonschema_validator():
    jsonschema = pytest.importorskip("jsonschema")
    schema = json.loads(published_schema_path().read_text(encoding="utf-8"))
    cls = jsonschema.validators.validator_for(schema)
    cls.check_schema(schema)
    return cls(schema, format_checker=cls.FORMAT_CHECKER)


@pytest.mark.parametrize("path", VALID, ids=[p.stem for p in VALID])
def test_valid_fixtures_pass_published_json_schema(path: Path) -> None:
    _jsonschema_validator().validate(_load(path))


@pytest.mark.parametrize("path", INVALID, ids=[p.stem for p in INVALID])
def test_invalid_fixtures_fail_published_json_schema(path: Path) -> None:
    validator = _jsonschema_validator()
    assert not validator.is_valid(_load(path))
