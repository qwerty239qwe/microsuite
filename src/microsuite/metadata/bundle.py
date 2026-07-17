"""Validation helpers for a referenced MicroSuite run bundle."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
from typing import Any

from microsuite.metadata.schemas import (
    READS_MANIFEST_VERSION,
    RESOLVED_CONFIG_VERSION,
    RUN_MANIFEST_VERSION,
    SCHEMA_VERSION,
    WORKFLOW_VERSION,
)
from microsuite.metadata.validate import validate_metadata


def sha256_file(path: Path) -> str:
    """Return the lowercase SHA-256 digest for one file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path, label: str) -> tuple[dict[str, Any] | None, list[str]]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return None, [f"{label}: cannot read JSON: {exc}"]
    if not isinstance(value, dict):
        return None, [f"{label}: expected a JSON object"]
    return value, []


def _resolve_internal(run_dir: Path, value: str, label: str) -> tuple[Path | None, list[str]]:
    posix = PurePosixPath(value)
    if "\\" in value or posix.is_absolute() or ".." in posix.parts or value in ("", "."):
        return None, [f"{label}: unsafe internal path {value!r}"]
    root = run_dir.resolve()
    target = root.joinpath(*posix.parts).resolve()
    if target != root and not target.is_relative_to(root):
        return None, [f"{label}: path escapes run bundle: {value!r}"]
    return target, []


def _validate_reference(
    run_dir: Path,
    reference: dict[str, Any],
    *,
    label: str,
    expected_schema: str,
) -> tuple[dict[str, Any] | None, list[str]]:
    errors: list[str] = []
    if reference.get("schema_version") != expected_schema:
        errors.append(
            f"{label}.schema_version: expected {expected_schema!r}, "
            f"got {reference.get('schema_version')!r}"
        )
    path, path_errors = _resolve_internal(run_dir, reference.get("path", ""), f"{label}.path")
    errors.extend(path_errors)
    if path is None:
        return None, errors
    if not path.is_file():
        errors.append(f"{label}: referenced file does not exist: {reference.get('path')!r}")
        return None, errors
    expected_digest = reference.get("sha256")
    if expected_digest is not None and sha256_file(path) != expected_digest:
        errors.append(f"{label}.sha256: checksum mismatch")
    document, load_errors = _load_json(path, label)
    errors.extend(load_errors)
    if document is not None:
        errors.extend(f"{label}: {error}" for error in validate_metadata(document, expected_schema))
    return document, errors


def validate_run_bundle(run_dir: Path) -> list[str]:
    """Validate a run manifest and every document it explicitly references."""
    run_dir = run_dir.resolve()
    manifest_path = run_dir / "run_manifest.json"
    manifest, errors = _load_json(manifest_path, "run_manifest.json")
    if manifest is None:
        return errors
    errors.extend(validate_metadata(manifest, RUN_MANIFEST_VERSION))
    if errors:
        return errors

    workflow, ref_errors = _validate_reference(
        run_dir,
        manifest["workflow"],
        label="workflow",
        expected_schema=WORKFLOW_VERSION,
    )
    errors.extend(ref_errors)
    _, ref_errors = _validate_reference(
        run_dir,
        manifest["resolved_config"],
        label="resolved_config",
        expected_schema=RESOLVED_CONFIG_VERSION,
    )
    errors.extend(ref_errors)
    if manifest.get("reads_manifest") is not None:
        _, ref_errors = _validate_reference(
            run_dir,
            manifest["reads_manifest"],
            label="reads_manifest",
            expected_schema=READS_MANIFEST_VERSION,
        )
        errors.extend(ref_errors)

    workflow_nodes = {node["id"] for node in workflow["nodes"]} if workflow else set()
    for i, stage in enumerate(manifest["stages"]):
        if workflow and stage["node_id"] not in workflow_nodes:
            errors.append(f"stages[{i}].node_id: absent from workflow: {stage['node_id']!r}")
        for j, reference in enumerate(stage["result_files"]):
            stage_result, ref_errors = _validate_reference(
                run_dir,
                reference,
                label=f"stages[{i}].result_files[{j}]",
                expected_schema=SCHEMA_VERSION,
            )
            errors.extend(ref_errors)
            if stage_result is not None and stage_result.get("run_id") != manifest["run_id"]:
                errors.append(
                    f"stages[{i}].result_files[{j}].run_id: does not match run manifest"
                )
    return errors
