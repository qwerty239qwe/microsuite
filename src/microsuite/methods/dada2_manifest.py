from __future__ import annotations

import json
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError

MANIFEST_FILENAME = "dada2_denoise_manifest.json"

_VERSION_KEYS = ("dada2_version", "r_version")


def read_r_params(path: Path) -> dict:
    if not path.exists():
        raise MicrobiomeSuiteError(f"DADA2 R params file was not written: {path}")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise MicrobiomeSuiteError(f"Could not parse DADA2 R params file {path}: {exc}") from exc


def build_manifest(r_params: dict, wrapper: dict) -> dict:
    params = {
        key: value
        for key, value in r_params.items()
        if key not in _VERSION_KEYS and value is not None
    }
    return {
        "tool": {key: r_params.get(key) for key in _VERSION_KEYS},
        "dada2_params": params,
        "run": dict(wrapper),
    }


def write_manifest(manifest: dict, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / MANIFEST_FILENAME
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return path
