from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class RefDbSource:
    name: str
    version: str
    locator: str | None = None


@dataclass(frozen=True)
class RefDbSpec:
    name: str
    version: str
    provider: str = "biodbs"
    target: str = "16S"
    build_targets: tuple[str, ...] = ("vsearch",)
    sources: tuple[RefDbSource, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class RawRefDb:
    sequences: Path
    taxonomy: Path
    qza: Path | None = None


@dataclass(frozen=True)
class BuiltArtifact:
    path: Path
    build_target: str
    checksum: str
