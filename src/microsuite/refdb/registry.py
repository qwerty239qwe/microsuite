from __future__ import annotations

import hashlib
import json
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.spec import BuiltArtifact


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


class RefDbRegistry:
    def __init__(self, root: Path) -> None:
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.root / "manifest.json"

    def _key(self, name: str, version: str, build_target: str) -> str:
        return f"{name}@{version}:{build_target}"

    def _load(self) -> dict[str, dict[str, str]]:
        if not self.manifest_path.exists():
            return {}
        try:
            return json.loads(self.manifest_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError) as exc:
            raise MicrobiomeSuiteError(
                f"Reference-DB manifest is corrupt: {self.manifest_path}. "
                "Delete it to rebuild the cache."
            ) from exc

    def resolve(self, name: str, version: str, build_target: str) -> BuiltArtifact | None:
        entry = self._load().get(self._key(name, version, build_target))
        if entry is None:
            return None
        path = Path(entry["path"])
        if not path.exists() or sha256_file(path) != entry["checksum"]:
            return None
        return BuiltArtifact(path=path, build_target=build_target, checksum=entry["checksum"])

    def record(
        self, name: str, version: str, artifact: BuiltArtifact, provider: str
    ) -> None:
        manifest = self._load()
        manifest[self._key(name, version, artifact.build_target)] = {
            "path": str(artifact.path),
            "checksum": artifact.checksum,
            "provider": provider,
        }
        self.manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
        )
