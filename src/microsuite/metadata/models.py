"""Typed declaration models for the ``stage-result.v1`` envelope.

A stage declares its inputs, outputs, and provenance files explicitly (never
guessed from argv or folder scans). Declared paths MUST be absolute so the writer
can decide, unambiguously, whether a path is inside ``run_dir`` (serialised
relative) or outside it (``external: true``) — it never has to guess whether a
relative path was based on the process CWD, ``run_command(cwd=...)``, or
``run_dir``.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


def _require_absolute(path: str | Path, what: str) -> None:
    if not Path(path).is_absolute():
        raise ValueError(f"{what} path must be absolute when declared, got relative: {path!r}")


@dataclass(frozen=True)
class ArtifactCount:
    """A declared count with an explicit unit (e.g. 1842 ``features``)."""

    value: int
    unit: str


@dataclass(frozen=True)
class Artifact:
    """A declared stage input or output file/directory.

    ``path`` must be absolute at declaration time. The writer fills ``external``
    (outside ``run_dir``), ``exists``, and ``bytes`` at serialisation.
    """

    label: str
    path: str | Path
    format: str | None = None
    kind: str | None = None
    count: ArtifactCount | None = None
    required: bool = True
    external: bool = False

    def __post_init__(self) -> None:
        _require_absolute(self.path, "Artifact")


@dataclass(frozen=True)
class ProvenanceFile:
    """A reference to a bespoke sidecar (dada2 manifest, ancombc provenance, …)."""

    kind: str
    path: str | Path
    required: bool = True

    def __post_init__(self) -> None:
        _require_absolute(self.path, "ProvenanceFile")


@dataclass(frozen=True)
class StageError:
    """Structured, already-redacted failure information for the envelope."""

    type: str
    message: str
