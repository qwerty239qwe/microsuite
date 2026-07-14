"""Dependency-free validator for the ``stage-result.v1`` envelope.

``validate`` checks structure against the internal schema grammar
(:mod:`microsuite.metadata.schemas`). ``validate_stage_result`` additionally
enforces the cross-field status/exit/alias/required invariants that are the
authoritative writer-side guard. Both return a list of human-readable error
strings; an empty list means valid.
"""

from __future__ import annotations

import math
from datetime import datetime
from typing import Any

from microsuite.metadata.schemas import SCHEMA_VERSION, SCHEMAS

_FAILING_SUB = {"failed", "launch_failed", "timed_out"}


def _is_rfc3339(value: str) -> bool:
    if not value.endswith("Z"):
        return False
    try:
        datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return False
    return True


def _check(value: Any, spec: dict[str, Any], path: str) -> list[str]:
    if value is None:
        if spec.get("nullable"):
            return []
        return [f"{path}: must not be null"]

    kind = spec.get("type")
    errors: list[str] = []

    if kind == "str":
        if not isinstance(value, str):
            return [f"{path}: expected str, got {type(value).__name__}"]
        if spec.get("format") == "rfc3339" and not _is_rfc3339(value):
            errors.append(f"{path}: not an RFC 3339 'Z' timestamp: {value!r}")
    elif kind == "int":
        if type(value) is not int:  # rejects bool
            return [f"{path}: expected int, got {type(value).__name__}"]
    elif kind == "number":
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return [f"{path}: expected number, got {type(value).__name__}"]
        if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
            return [f"{path}: NaN/inf not allowed"]
    elif kind == "bool":
        if not isinstance(value, bool):
            return [f"{path}: expected bool, got {type(value).__name__}"]
    elif kind == "object":
        if not isinstance(value, dict):
            return [f"{path}: expected object, got {type(value).__name__}"]
        return _check_object(value, spec, path)
    elif kind == "array":
        if not isinstance(value, list):
            return [f"{path}: expected array, got {type(value).__name__}"]
        item_spec = spec.get("items")
        if item_spec:
            for i, item in enumerate(value):
                errors.extend(_check(item, item_spec, f"{path}[{i}]"))
        return errors
    else:  # pragma: no cover - schema authoring error
        return [f"{path}: unknown schema type {kind!r}"]

    if "const" in spec and value != spec["const"]:
        errors.append(f"{path}: expected const {spec['const']!r}, got {value!r}")
    if "enum" in spec and value not in spec["enum"]:
        errors.append(f"{path}: {value!r} not in {spec['enum']}")
    if "min" in spec and isinstance(value, (int, float)) and value < spec["min"]:
        errors.append(f"{path}: {value} < min {spec['min']}")
    return errors


def _check_object(value: dict[str, Any], spec: dict[str, Any], path: str) -> list[str]:
    errors: list[str] = []
    fields = spec.get("fields", {})
    for key in spec.get("required", []):
        if key not in value:
            errors.append(f"{path}: missing required field {key!r}")
    if not spec.get("allow_unknown", False):
        for key in value:
            if key not in fields:
                errors.append(f"{path}: unknown field {key!r}")
    for key, sub_spec in fields.items():
        if key in value:
            sub_path = f"{path}.{key}" if path else key
            errors.extend(_check(value[key], sub_spec, sub_path))
    return errors


def validate(payload: Any, schema_name: str = SCHEMA_VERSION) -> list[str]:
    """Structural validation against the named internal schema."""
    schema = SCHEMAS.get(schema_name)
    if schema is None:
        return [f"unknown schema {schema_name!r}"]
    if not isinstance(payload, dict):
        return [f"payload: expected object, got {type(payload).__name__}"]
    return _check_object(payload, schema, "")


def expected_alias(payload: dict[str, Any]) -> tuple[Any, Any]:
    """The (command, exit_code) the top-level alias must equal for this status."""
    status = payload.get("status")
    subs = payload.get("subprocesses") or []
    if status == "completed":
        chosen = subs[-1] if subs else None
    elif status in ("failed", "timed_out"):
        responsible = [s for s in subs if s.get("status") in _FAILING_SUB]
        chosen = responsible[-1] if responsible else None
    else:  # cancelled / running
        chosen = None
    if chosen is None:
        return None, None
    return chosen.get("command"), chosen.get("exit_code")


def _invariants(payload: dict[str, Any]) -> list[str]:
    errors: list[str] = []
    status = payload["status"]
    exit_code = payload.get("exit_code")
    error = payload.get("error")
    subs = payload.get("subprocesses") or []

    if status == "completed":
        if error is not None:
            errors.append("completed: error must be null")
        if exit_code not in (0, None):
            errors.append("completed: exit_code must be 0 or null")
        for i, sp in enumerate(subs):
            if sp.get("required", True) and sp.get("status") != "completed":
                errors.append(f"completed: required subprocess[{i}] did not complete")
        for group in ("outputs", "provenance_files"):
            for i, art in enumerate(payload.get(group, [])):
                if art.get("required", True) and art.get("exists") is not True:
                    errors.append(f"completed: required {group}[{i}] does not exist")
    elif status == "failed":
        if error is None:
            errors.append("failed: error must be non-null")
    elif status == "timed_out":
        if exit_code is not None:
            errors.append("timed_out: exit_code must be null")
        if error is None:
            errors.append("timed_out: error must be non-null")
        if not any(s.get("status") == "timed_out" for s in subs):
            errors.append("timed_out: needs a timed_out subprocess")
    elif status == "cancelled":
        if exit_code is not None:
            errors.append("cancelled: exit_code must be null")
        if error is None:
            errors.append("cancelled: error must be non-null")

    # Per-subprocess status/exit consistency.
    for i, sp in enumerate(subs):
        sp_status, sp_exit = sp.get("status"), sp.get("exit_code")
        if sp_status == "completed" and sp_exit != 0:
            errors.append(f"subprocess[{i}] completed requires exit_code 0")
        if sp_status == "failed" and not (isinstance(sp_exit, int) and sp_exit != 0):
            errors.append(f"subprocess[{i}] failed requires a non-zero int exit_code")
        if sp_status in ("timed_out", "launch_failed") and sp_exit is not None:
            errors.append(f"subprocess[{i}] {sp_status} requires null exit_code")

    # Alias consistency.
    exp_cmd, exp_exit = expected_alias(payload)
    if payload.get("command") != exp_cmd:
        errors.append("command alias does not match the responsible subprocess")
    if payload.get("exit_code") != exp_exit:
        errors.append("exit_code alias does not match the responsible subprocess")
    return errors


def validate_stage_result(payload: Any) -> list[str]:
    """Structural + semantic-invariant validation of a stage-result payload."""
    errors = validate(payload, SCHEMA_VERSION)
    if errors:
        return errors
    return _invariants(payload)
