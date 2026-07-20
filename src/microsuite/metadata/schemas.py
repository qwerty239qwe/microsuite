"""Internal schema definitions for versioned MicroSuite metadata documents.

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

RESOLVED_CONFIG_VERSION = "resolved-config.v1"
WORKFLOW_VERSION = "workflow.v1"
RUN_MANIFEST_VERSION = "run-manifest.v1"
READS_MANIFEST_VERSION = "reads-manifest.v1"

RESOLVED_CONFIG_V1 = {
    "type": "object",
    "allow_unknown": True,
    "required": ["schema_version", "generated_at", "producer", "config"],
    "fields": {
        "schema_version": {"type": "str", "const": RESOLVED_CONFIG_VERSION},
        "generated_at": {"type": "str", "format": "rfc3339"},
        "producer": {
            "type": "object",
            "required": ["name", "version"],
            "allow_unknown": True,
            "fields": {"name": {"type": "str"}, "version": {"type": "str"}},
        },
        "config": _OBJECT_ANY,
    },
}

_PRODUCER = {
    "type": "object",
    "required": ["name", "version"],
    "allow_unknown": True,
    "fields": {"name": {"type": "str"}, "version": {"type": "str"}},
}

_DOCUMENT_REF = {
    "type": "object",
    "required": ["path", "schema_version"],
    "allow_unknown": True,
    "fields": {
        "path": {"type": "str"},
        "schema_version": {"type": "str"},
        "sha256": {"type": "str", "nullable": True},
    },
}

_WORKFLOW_NODE = {
    "type": "object",
    "required": ["id", "stage_type", "process_name"],
    "allow_unknown": True,
    "fields": {
        "id": {"type": "str"},
        "stage_type": {"type": "str"},
        "process_name": {"type": "str"},
        "method": {"type": "str", "nullable": True},
        "branch": {"type": "str", "nullable": True},
        "rank": {"type": "str", "nullable": True},
        "inputs": {"type": "array", "items": {"type": "str"}},
        "outputs": {"type": "array", "items": {"type": "str"}},
    },
}

_WORKFLOW_EDGE = {
    "type": "object",
    "required": ["source", "target"],
    "allow_unknown": True,
    "fields": {
        "source": {"type": "str"},
        "target": {"type": "str"},
        "kind": {"type": "str", "nullable": True},
    },
}

_WORKFLOW_BRANCH = {
    "type": "object",
    "required": ["id", "enabled", "node_ids"],
    "allow_unknown": True,
    "fields": {
        "id": {"type": "str"},
        "enabled": {"type": "bool"},
        "node_ids": {"type": "array", "items": {"type": "str"}},
    },
}

_WORKFLOW_METHOD = {
    "type": "object",
    "required": ["id", "name"],
    "allow_unknown": True,
    "fields": {
        "id": {"type": "str"},
        "name": {"type": "str"},
        "version": {"type": "str", "nullable": True},
        "branch": {"type": "str", "nullable": True},
    },
}

_INTERMEDIATE_FILE = {
    "type": "object",
    "required": ["id", "label", "kind", "format"],
    "allow_unknown": True,
    "fields": {
        "id": {"type": "str"},
        "label": {"type": "str"},
        "kind": {"type": "str"},
        "format": {"type": "str"},
        "path": {"type": "str", "nullable": True},
        "producer_node": {"type": "str", "nullable": True},
        "consumer_nodes": {"type": "array", "items": {"type": "str"}},
        "retained": {"type": "bool"},
    },
}

WORKFLOW_V1 = {
    "type": "object",
    "allow_unknown": True,
    "required": [
        "schema_version",
        "workflow_id",
        "workflow_version",
        "workflow_hash",
        "generated_at",
        "producer",
        "engine",
        "nodes",
        "edges",
        "branches",
        "methods",
        "intermediate_files",
    ],
    "fields": {
        "schema_version": {"type": "str", "const": WORKFLOW_VERSION},
        "workflow_id": {"type": "str"},
        "workflow_version": {"type": "str"},
        "workflow_hash": {"type": "str"},
        "generated_at": {"type": "str", "format": "rfc3339"},
        "producer": _PRODUCER,
        "engine": {
            "type": "object",
            "required": ["name", "version"],
            "allow_unknown": True,
            "fields": {"name": {"type": "str"}, "version": {"type": "str"}},
        },
        "nodes": {"type": "array", "items": _WORKFLOW_NODE},
        "edges": {"type": "array", "items": _WORKFLOW_EDGE},
        "branches": {"type": "array", "items": _WORKFLOW_BRANCH},
        "methods": {"type": "array", "items": _WORKFLOW_METHOD},
        "intermediate_files": {"type": "array", "items": _INTERMEDIATE_FILE},
        "native_dag": {"type": "object", "nullable": True, "allow_unknown": True, "fields": {}},
    },
}

_STAGE_SUMMARY = {
    "type": "object",
    "required": [
        "node_id",
        "status",
        "expected_tasks",
        "completed_tasks",
        "failed_tasks",
        "result_files",
    ],
    "allow_unknown": True,
    "fields": {
        "node_id": {"type": "str"},
        "status": {"type": "str", "enum": _STATUS_ENUM + ["pending", "skipped"]},
        "expected_tasks": {"type": "int", "min": 0},
        "completed_tasks": {"type": "int", "min": 0},
        "failed_tasks": {"type": "int", "min": 0},
        "result_files": {"type": "array", "items": _DOCUMENT_REF},
    },
}

_PROGRESS = {
    "type": "object",
    "required": ["total", "completed", "failed", "running", "pending", "fraction"],
    "allow_unknown": True,
    "fields": {
        "total": {"type": "int", "min": 0},
        "completed": {"type": "int", "min": 0},
        "failed": {"type": "int", "min": 0},
        "running": {"type": "int", "min": 0},
        "pending": {"type": "int", "min": 0},
        "fraction": {"type": "number", "min": 0, "max": 1},
    },
}

RUN_MANIFEST_V1 = {
    "type": "object",
    "allow_unknown": True,
    "required": [
        "schema_version",
        "run_id",
        "dataset_id",
        "analysis_id",
        "status",
        "created_at",
        "updated_at",
        "producer",
        "workflow",
        "resolved_config",
        "progress",
        "stages",
        "artifacts",
        "summary_metrics",
        "provenance",
    ],
    "fields": {
        "schema_version": {"type": "str", "const": RUN_MANIFEST_VERSION},
        "run_id": {"type": "str"},
        "dataset_id": {"type": "str"},
        "analysis_id": {"type": "str"},
        "status": {
            "type": "str",
            "enum": ["initializing", "running", "completed", "failed", "cancelled", "timed_out"],
        },
        "created_at": {"type": "str", "format": "rfc3339"},
        "updated_at": {"type": "str", "format": "rfc3339"},
        "producer": _PRODUCER,
        "workflow": _DOCUMENT_REF,
        "resolved_config": _DOCUMENT_REF,
        "reads_manifest": {
            "type": "object",
            "nullable": True,
            "allow_unknown": True,
            "fields": _DOCUMENT_REF["fields"],
            "required": _DOCUMENT_REF["required"],
        },
        "progress": _PROGRESS,
        "stages": {"type": "array", "items": _STAGE_SUMMARY},
        "artifacts": {"type": "array", "items": _ARTIFACT},
        "summary_metrics": _OBJECT_ANY,
        "provenance": _OBJECT_ANY,
        "error": {"type": "object", "nullable": True, "allow_unknown": True, "fields": {}},
    },
}

_READ_FILE = {
    "type": "object",
    "required": ["path", "read", "external"],
    "allow_unknown": True,
    "fields": {
        "path": {"type": "str"},
        "read": {"type": "str", "enum": ["R1", "R2", "single"]},
        "external": {"type": "bool"},
        "bytes": {"type": "int", "min": 0, "nullable": True},
        "sha256": {"type": "str", "nullable": True},
    },
}

_READ_SAMPLE = {
    "type": "object",
    "required": ["sample_id", "layout", "files"],
    "allow_unknown": True,
    "fields": {
        "sample_id": {"type": "str"},
        "layout": {"type": "str", "enum": ["SE", "PE"]},
        "files": {"type": "array", "items": _READ_FILE},
        "metadata": _OBJECT_ANY,
    },
}

READS_MANIFEST_V1 = {
    "type": "object",
    "allow_unknown": True,
    "required": [
        "schema_version",
        "dataset_id",
        "generated_at",
        "producer",
        "layout",
        "samples",
    ],
    "fields": {
        "schema_version": {"type": "str", "const": READS_MANIFEST_VERSION},
        "dataset_id": {"type": "str"},
        "generated_at": {"type": "str", "format": "rfc3339"},
        "producer": _PRODUCER,
        "layout": {"type": "str", "enum": ["SE", "PE", "mixed"]},
        "samples": {"type": "array", "items": _READ_SAMPLE},
        "selection": _OBJECT_ANY,
        "source": _OBJECT_ANY,
    },
}

SCHEMAS = {
    SCHEMA_VERSION: STAGE_RESULT_V1,
    RESOLVED_CONFIG_VERSION: RESOLVED_CONFIG_V1,
    WORKFLOW_VERSION: WORKFLOW_V1,
    RUN_MANIFEST_VERSION: RUN_MANIFEST_V1,
    READS_MANIFEST_VERSION: READS_MANIFEST_V1,
}


def published_schema_path(schema_name: str = SCHEMA_VERSION) -> Path:
    """Path to a packaged, language-neutral JSON Schema (draft 2020-12)."""
    return Path(str(files("microsuite.metadata._schema").joinpath(f"{schema_name}.schema.json")))
