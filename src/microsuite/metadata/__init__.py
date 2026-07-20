"""Structured run metadata for Microboard (sub-project A).

Public surface: typed declaration models, the ``stage_execution`` boundary that
finalizes exactly one ``stage-result.v1`` envelope per stage, and the schema
registry + validator that guard the envelope contract. See
``docs/superpowers/specs/2026-07-14-stage-result-metadata-A-design.md``.
"""

from __future__ import annotations

from microsuite.metadata.bundle import sha256_file, validate_run_bundle
from microsuite.metadata.config import write_resolved_config
from microsuite.metadata.context import WorkflowContext, workflow_context_from_env
from microsuite.metadata.documents import (
    write_metadata_document,
    write_reads_manifest,
    write_run_manifest,
    write_workflow,
)
from microsuite.metadata.models import (
    Artifact,
    ArtifactCount,
    ProvenanceFile,
    StageError,
)
from microsuite.metadata.stage import (
    StageRecord,
    active_stage,
    stage_execution,
)
from microsuite.metadata.validate import (
    validate,
    validate_metadata,
    validate_reads_manifest,
    validate_run_manifest,
    validate_stage_result,
    validate_workflow,
)

__all__ = [
    "Artifact",
    "ArtifactCount",
    "ProvenanceFile",
    "StageError",
    "StageRecord",
    "WorkflowContext",
    "active_stage",
    "sha256_file",
    "stage_execution",
    "validate",
    "validate_metadata",
    "validate_reads_manifest",
    "validate_run_bundle",
    "validate_run_manifest",
    "validate_stage_result",
    "validate_workflow",
    "workflow_context_from_env",
    "write_metadata_document",
    "write_reads_manifest",
    "write_resolved_config",
    "write_run_manifest",
    "write_workflow",
]
