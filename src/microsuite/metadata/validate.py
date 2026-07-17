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
from pathlib import PurePosixPath
from typing import Any

from microsuite.metadata.schemas import (
    READS_MANIFEST_VERSION,
    RESOLVED_CONFIG_VERSION,
    RUN_MANIFEST_VERSION,
    SCHEMA_VERSION,
    SCHEMAS,
    WORKFLOW_VERSION,
)

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
    if "max" in spec and isinstance(value, (int, float)) and value > spec["max"]:
        errors.append(f"{path}: {value} > max {spec['max']}")
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


def _duplicate_values(values: list[Any]) -> set[Any]:
    seen: set[Any] = set()
    duplicates: set[Any] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    return duplicates


def _validate_bundle_path(value: str, path: str) -> list[str]:
    if "\\" in value:
        return [f"{path}: internal path must use '/' separators"]
    candidate = PurePosixPath(value)
    if candidate.is_absolute() or ".." in candidate.parts or value in ("", "."):
        return [f"{path}: internal path must be a safe relative path"]
    return []


def validate_workflow(payload: Any) -> list[str]:
    """Validate workflow structure and graph-reference invariants."""
    errors = validate(payload, WORKFLOW_VERSION)
    if errors:
        return errors
    node_ids = [node["id"] for node in payload["nodes"]]
    node_set = set(node_ids)
    for duplicate in sorted(_duplicate_values(node_ids)):
        errors.append(f"nodes: duplicate id {duplicate!r}")
    for i, edge in enumerate(payload["edges"]):
        for endpoint in ("source", "target"):
            if edge[endpoint] not in node_set:
                errors.append(f"edges[{i}].{endpoint}: unknown node {edge[endpoint]!r}")
    for i, branch in enumerate(payload["branches"]):
        for node_id in branch["node_ids"]:
            if node_id not in node_set:
                errors.append(f"branches[{i}].node_ids: unknown node {node_id!r}")
    for i, file_def in enumerate(payload["intermediate_files"]):
        producer = file_def.get("producer_node")
        if producer is not None and producer not in node_set:
            errors.append(f"intermediate_files[{i}].producer_node: unknown node {producer!r}")
        for consumer in file_def.get("consumer_nodes", []):
            if consumer not in node_set:
                errors.append(
                    f"intermediate_files[{i}].consumer_nodes: unknown node {consumer!r}"
                )
    return errors


def validate_reads_manifest(payload: Any) -> list[str]:
    """Validate normalized sample/read layout invariants."""
    errors = validate(payload, READS_MANIFEST_VERSION)
    if errors:
        return errors
    sample_ids = [sample["sample_id"] for sample in payload["samples"]]
    for duplicate in sorted(_duplicate_values(sample_ids)):
        errors.append(f"samples: duplicate sample_id {duplicate!r}")
    layouts: set[str] = set()
    for i, sample in enumerate(payload["samples"]):
        layout = sample["layout"]
        layouts.add(layout)
        reads = [item["read"] for item in sample["files"]]
        if layout == "PE" and sorted(reads) != ["R1", "R2"]:
            errors.append(f"samples[{i}]: PE layout requires exactly one R1 and one R2")
        if layout == "SE" and len(reads) != 1:
            errors.append(f"samples[{i}]: SE layout requires exactly one read file")
        if layout == "SE" and reads and reads[0] not in ("R1", "single"):
            errors.append(f"samples[{i}]: SE read must be labelled 'single' or 'R1'")
    declared = payload["layout"]
    expected = next(iter(layouts)) if len(layouts) == 1 else "mixed"
    if layouts and declared != expected:
        errors.append(f"layout: declared {declared!r}, expected {expected!r} from samples")
    return errors


def validate_run_manifest(payload: Any) -> list[str]:
    """Validate run progress and internal reference invariants."""
    errors = validate(payload, RUN_MANIFEST_VERSION)
    if errors:
        return errors
    for name in ("workflow", "resolved_config", "reads_manifest"):
        reference = payload.get(name)
        if reference is not None:
            errors.extend(_validate_bundle_path(reference["path"], f"{name}.path"))
    stage_ids = [stage["node_id"] for stage in payload["stages"]]
    for duplicate in sorted(_duplicate_values(stage_ids)):
        errors.append(f"stages: duplicate node_id {duplicate!r}")
    for i, stage in enumerate(payload["stages"]):
        for j, reference in enumerate(stage["result_files"]):
            errors.extend(
                _validate_bundle_path(reference["path"], f"stages[{i}].result_files[{j}].path")
            )
    progress = payload["progress"]
    accounted = sum(progress[key] for key in ("completed", "failed", "running", "pending"))
    if accounted != progress["total"]:
        errors.append(
            f"progress: task counts sum to {accounted}, expected total {progress['total']}"
        )
    expected_fraction = (
        1.0
        if progress["total"] == 0 and payload["status"] == "completed"
        else 0.0
        if progress["total"] == 0
        else (progress["completed"] + progress["failed"]) / progress["total"]
    )
    if not math.isclose(progress["fraction"], expected_fraction, rel_tol=0, abs_tol=1e-9):
        errors.append(
            f"progress.fraction: {progress['fraction']} does not match terminal/total "
            f"fraction {expected_fraction}"
        )
    return errors


def validate_metadata(payload: Any, schema_name: str | None = None) -> list[str]:
    """Validate one supported metadata document, including cross-field rules."""
    if schema_name is None:
        if not isinstance(payload, dict):
            return [f"payload: expected object, got {type(payload).__name__}"]
        schema_name = payload.get("schema_version")
    if schema_name == SCHEMA_VERSION:
        return validate_stage_result(payload)
    if schema_name == WORKFLOW_VERSION:
        return validate_workflow(payload)
    if schema_name == RUN_MANIFEST_VERSION:
        return validate_run_manifest(payload)
    if schema_name == READS_MANIFEST_VERSION:
        return validate_reads_manifest(payload)
    if schema_name == RESOLVED_CONFIG_VERSION:
        return validate(payload, RESOLVED_CONFIG_VERSION)
    return [f"unknown schema {schema_name!r}"]
