"""Atomic writers for versioned workflow-level metadata documents."""

from __future__ import annotations

import time
from collections.abc import Mapping
from pathlib import Path
from typing import Any

from microsuite import __version__
from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata.redact import json_safe, redact_params
from microsuite.metadata.schemas import (
    READS_MANIFEST_VERSION,
    RUN_MANIFEST_VERSION,
    WORKFLOW_VERSION,
)
from microsuite.metadata.stage import _atomic_write
from microsuite.metadata.validate import validate_metadata


def _utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def write_metadata_document(
    target: Path,
    payload: Mapping[str, Any],
    *,
    schema_name: str,
    redact: bool = True,
) -> Path:
    """Validate and atomically write one metadata document.

    The caller supplies the semantic payload. This function enforces the schema
    version, converts values to JSON-safe forms, redacts secret-bearing fields,
    and refuses to publish invalid data.
    """
    document = json_safe(dict(payload))
    document["schema_version"] = schema_name
    if redact:
        document, _ = redact_params(document)
    errors = validate_metadata(document, schema_name)
    if errors:
        raise MicrobiomeSuiteError(f"Refusing to write invalid {schema_name} metadata: {errors}")
    target.parent.mkdir(parents=True, exist_ok=True)
    _atomic_write(target, document)
    return target


def write_workflow(target: Path, payload: Mapping[str, Any]) -> Path:
    """Write a ``workflow.v1`` snapshot, filling producer/timestamp defaults."""
    document = dict(payload)
    document.setdefault("generated_at", _utc_now())
    document.setdefault("producer", {"name": "microsuite", "version": __version__})
    return write_metadata_document(target, document, schema_name=WORKFLOW_VERSION)


def write_reads_manifest(target: Path, payload: Mapping[str, Any]) -> Path:
    """Write a normalized ``reads-manifest.v1`` document."""
    document = dict(payload)
    document.setdefault("generated_at", _utc_now())
    document.setdefault("producer", {"name": "microsuite", "version": __version__})
    return write_metadata_document(target, document, schema_name=READS_MANIFEST_VERSION)


def write_run_manifest(target: Path, payload: Mapping[str, Any]) -> Path:
    """Write a ``run-manifest.v1`` document."""
    document = dict(payload)
    document.setdefault("producer", {"name": "microsuite", "version": __version__})
    return write_metadata_document(target, document, schema_name=RUN_MANIFEST_VERSION)
