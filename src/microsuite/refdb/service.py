# src/microsuite/refdb/service.py
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite.refdb.build import merge_raw
from microsuite.refdb.paths import refdb_cache_dir
from microsuite.refdb.providers import get_provider
from microsuite.refdb.registry import RefDbRegistry
from microsuite.refdb.spec import BuiltArtifact, RefDbSpec


def _work_dir(registry: RefDbRegistry, spec: RefDbSpec, build_target: str) -> Path:
    work = registry.root / f"{spec.name}@{spec.version}" / build_target
    work.mkdir(parents=True, exist_ok=True)
    return work


def fetch_refdb(
    spec: RefDbSpec,
    build_target: str,
    *,
    force: bool = False,
    registry: RefDbRegistry | None = None,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> BuiltArtifact:
    registry = registry or RefDbRegistry(refdb_cache_dir())
    if not force:
        cached = registry.resolve(spec.name, spec.version, build_target)
        if cached is not None:
            return cached

    provider = get_provider(spec.provider)
    work = _work_dir(registry, spec, build_target)
    if spec.sources:
        raws = []
        for source in spec.sources:
            sub = replace(spec, name=source.name, version=source.version, sources=())
            sub_dir = work / "sources" / source.name
            sub_dir.mkdir(parents=True, exist_ok=True)
            raws.append(provider.fetch(sub, out_dir=sub_dir))
        raw = merge_raw(raws, out_dir=work / "merged") if len(raws) > 1 else raws[0]
    else:
        raw = provider.fetch(spec, out_dir=work / "fetch")

    artifact = provider.build(
        raw, build_target, out_dir=work / "build", run_dir=run_dir, timeout=timeout
    )
    registry.record(spec.name, spec.version, artifact, spec.provider)
    return artifact


def resolve_classifier(value: str, *, registry: RefDbRegistry | None = None) -> Path:
    if not value.startswith("refdb:"):
        return Path(value)
    body = value[len("refdb:") :]
    build_target = "vsearch"
    if ":" in body:
        body, build_target = body.rsplit(":", 1)
    if "@" not in body:
        raise MicrobiomeSuiteError(
            f"Malformed refdb reference '{value}'. Expected refdb:<name>@<version>[:<build>]."
        )
    name, version = body.split("@", 1)
    registry = registry or RefDbRegistry(refdb_cache_dir())
    art = registry.resolve(name, version, build_target)
    if art is None:
        raise MicrobiomeSuiteError(
            f"Reference DB '{name}@{version}:{build_target}' is not in the cache. "
            f"Run: microsuite refdb fetch {name} --version {version} --build {build_target}"
        )
    return art.path
