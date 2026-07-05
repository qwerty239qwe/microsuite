from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from microsuite.refdb.build import build_artifact
from microsuite.refdb.spec import BuiltArtifact, RawRefDb, RefDbSpec


class RefDbProvider(ABC):
    name: str

    @abstractmethod
    def fetch(self, spec: RefDbSpec, out_dir: Path) -> RawRefDb: ...

    def build(
        self,
        raw: RawRefDb,
        build_target: str,
        out_dir: Path,
        *,
        force: bool = False,
        run_dir: Path | None = None,
        timeout: float | None = None,
    ) -> BuiltArtifact:
        return build_artifact(
            raw, build_target, out_dir, force=force, run_dir=run_dir, timeout=timeout
        )
