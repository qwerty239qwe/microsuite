"""Internal schema definitions for the ``stage-result.v1`` envelope.

This module holds the *MicroSuite-specific* schema format consumed by the
dependency-free Python validator in :mod:`microsuite.metadata.validate`. It is
**not** the published, language-neutral interoperability contract — that is the
draft 2020-12 JSON Schema published at
``metadata/schemas/stage-result.v1.schema.json`` for Microboard/Nextflow/other
consumers. The two are kept in parity by tests.

Schema grammar (a small subset inspired by JSON Schema):
``{"type": str|int|number|bool|object|array}`` with optional ``nullable``,
``const``, ``enum``, ``min`` (numeric), ``format`` ("rfc3339"); objects carry
``required`` (list), ``allow_unknown`` (bool), ``fields`` (name -> spec); arrays
carry ``items`` (spec).
"""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path

SCHEMA_VERSION = "stage-result.v1"
_STATUS_ENUM = ["running", "completed", "failed", "timed_out", "cancelled"]
_SUBPROC_STATUS_ENUM = ["completed", "failed", "timed_out", "launch_failed"]

_COUNT = {
    "type": "object",
    "nullable": True,
    "required": ["value", "unit"],
    "allow_unknown": True,
    "fields": {"value": {"type": "int", "min": 0}, "unit": {"type": "str"}},
}

_ARTIFACT = {
    "type": "object",
    "required": ["label", "path"],
    "allow_unknown": True,
    "fields": {
        "label": {"type": "str"},
        "path": {"type": "str"},
        "format": {"type": "str", "nullable": True},
        "kind": {"type": "str", "nullable": True},
        "count": _COUNT,
        "required": {"type": "bool"},
        "external": {"type": "bool"},
        "exists": {"type": "bool"},
        "bytes": {"type": "int", "min": 0, "nullable": True},
    },
}

_PROVENANCE = {
    "type": "object",
    "required": ["kind", "path"],
    "allow_unknown": True,
    "fields": {
        "kind": {"type": "str"},
        "path": {"type": "str"},
        "required": {"type": "bool"},
        "external": {"type": "bool"},
        "exists": {"type": "bool"},
    },
}

_SUBPROCESS = {
    "type": "object",
    "required": ["command", "status", "exit_code", "duration_sec"],
    "allow_unknown": True,
    "fields": {
        "command": {"type": "array", "items": {"type": "str"}},
        "status": {"type": "str", "enum": _SUBPROC_STATUS_ENUM},
        "exit_code": {"type": "int", "nullable": True},
        "duration_sec": {"type": "number", "min": 0},
        "required": {"type": "bool"},
    },
}

_TIMING = {
    "type": "object",
    "required": ["started_at", "finished_at", "duration_sec"],
    "allow_unknown": True,
    "fields": {
        "started_at": {"type": "str", "format": "rfc3339"},
        "finished_at": {"type": "str", "format": "rfc3339"},
        "duration_sec": {"type": "number", "min": 0},
    },
}

_OBJECT_ANY = {"type": "object", "allow_unknown": True, "fields": {}}

STAGE_RESULT_V1 = {
    "type": "object",
    "allow_unknown": True,
    "required": [
        "schema_version",
        "run_id",
        "stage_run_id",
        "attempt",
        "stage",
        "task",
        "backend",
        "status",
        "exit_code",
        "error",
        "timing",
        "command",
        "subprocesses",
        "params",
        "inputs",
        "outputs",
        "provenance_files",
        "metrics",
        "software",
        "reference_db",
        "producer",
    ],
    "fields": {
        "schema_version": {"type": "str", "const": SCHEMA_VERSION},
        "run_id": {"type": "str"},
        "stage_run_id": {"type": "str"},
        "attempt": {"type": "int", "min": 1},
        "stage": {"type": "str"},
        "task": {"type": "str"},
        "backend": {"type": "str", "nullable": True},
        "status": {"type": "str", "enum": _STATUS_ENUM},
        "exit_code": {"type": "int", "nullable": True},
        "error": {
            "type": "object",
            "nullable": True,
            "required": ["type", "message"],
            "allow_unknown": True,
            "fields": {"type": {"type": "str"}, "message": {"type": "str"}},
        },
        "timing": _TIMING,
        "command": {"type": "array", "nullable": True, "items": {"type": "str"}},
        "subprocesses": {"type": "array", "items": _SUBPROCESS},
        "params": _OBJECT_ANY,
        "inputs": {"type": "array", "items": _ARTIFACT},
        "outputs": {"type": "array", "items": _ARTIFACT},
        "provenance_files": {"type": "array", "items": _PROVENANCE},
        "metrics": _OBJECT_ANY,
        "software": _OBJECT_ANY,
        "reference_db": {"type": "object", "nullable": True, "allow_unknown": True, "fields": {}},
        "producer": {
            "type": "object",
            "required": ["name", "version"],
            "allow_unknown": True,
            "fields": {"name": {"type": "str"}, "version": {"type": "str"}},
        },
        "workflow_id": {"type": "str", "nullable": True},
        "workflow_run_id": {"type": "str", "nullable": True},
        "dataset_id": {"type": "str", "nullable": True},
    },
}

SCHEMAS = {SCHEMA_VERSION: STAGE_RESULT_V1}


def published_schema_path() -> Path:
    """Path to the packaged, language-neutral JSON Schema (draft 2020-12)."""
    return Path(str(files("microsuite.metadata._schema").joinpath("stage-result.v1.schema.json")))
