from __future__ import annotations

from abc import ABC, abstractmethod
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.build import build_artifact
from microsuite.refdb.registry import sha256_file
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
        if raw.qza is not None:
            if build_target == "qiime2":
                return BuiltArtifact(raw.qza, "qiime2", sha256_file(raw.qza))
            if build_target in ("vsearch", "blast"):
                raise MicrobiomeSuiteError(
                    f"Provider '{getattr(self, 'name', type(self).__name__)}' returned a "
                    "pre-packaged QIIME2 '.qza' artifact, which cannot be converted into a "
                    f"'{build_target}' reference here. Use '--build qiime2' for this "
                    "provider, or choose a provider that yields raw FASTA/taxonomy files "
                    "if you need a vsearch or blast reference."
                )
        return build_artifact(
            raw, build_target, out_dir, force=force, run_dir=run_dir, timeout=timeout
        )
