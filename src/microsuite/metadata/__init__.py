from __future__ import annotations

"""Structured run metadata for Microboard (sub-project A).

Public surface: typed declaration models, the ``stage_execution`` boundary that
finalizes exactly one ``stage-result.v1`` envelope per stage, and the schema
registry + validator that guard the envelope contract. See
``docs/superpowers/specs/2026-07-14-stage-result-metadata-A-design.md``.
"""

from microsuite.metadata.models import (
    Artifact,
    ArtifactCount,
    ProvenanceFile,
    StageError,
)

__all__ = [
    "Artifact",
    "ArtifactCount",
    "ProvenanceFile",
    "StageError",
]
