from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

from microsuite import __version__

RESULTS_MANIFEST = "microsuite-results.json"
SCHEMA_VERSION = "microsuite-results.v1"


def write_results_manifest(
    run_dir: Path,
    *,
    command: list[str],
    log: dict[str, Any],
    started_at: float,
    duration_sec: float | None = None,
    exit_code: int | None = None,
) -> None:
    manifest_path = run_dir / RESULTS_MANIFEST
    manifest = _read_manifest(manifest_path, run_dir)
    execution = {
        "task": log.get("task"),
        "backend": log.get("backend"),
        "command": command,
        "started_at": started_at,
        "inputs": log.get("inputs"),
        "outputs": log.get("outputs"),
        "params": log.get("params"),
    }
    if duration_sec is not None:
        execution["duration_sec"] = duration_sec
    if exit_code is not None:
        execution["exit_code"] = exit_code
    execution = {key: value for key, value in execution.items() if value is not None}

    manifest["executions"].append(execution)
    manifest["artifacts"] = _merge_artifacts(
        manifest["artifacts"],
        _artifacts_from_outputs(log),
    )
    manifest["updated_at"] = _utc_now()
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )


def _read_manifest(path: Path, run_dir: Path) -> dict[str, Any]:
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    return {
        "schema_version": SCHEMA_VERSION,
        "producer": {"name": "microsuite", "version": __version__},
        "run_id": run_dir.name,
        "run_dir": str(run_dir),
        "created_at": _utc_now(),
        "updated_at": _utc_now(),
        "executions": [],
        "artifacts": [],
    }


def _artifacts_from_outputs(log: dict[str, Any]) -> list[dict[str, Any]]:
    outputs = log.get("outputs")
    if not isinstance(outputs, dict):
        return []
    artifacts: list[dict[str, Any]] = []
    for label, value in outputs.items():
        if value is None:
            continue
        path = str(value)
        artifacts.append(
            {
                "id": _artifact_id(log.get("task"), log.get("backend"), label, path),
                "kind": _artifact_kind(log.get("task"), label),
                "label": label,
                "path": path,
                "format": _format_from_path(path),
                "task": log.get("task"),
                "backend": log.get("backend"),
            }
        )
    return artifacts


def _merge_artifacts(
    existing: list[dict[str, Any]], new: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    merged = {artifact["id"]: artifact for artifact in existing}
    for artifact in new:
        merged[artifact["id"]] = artifact
    return list(merged.values())


def _artifact_id(task: object, backend: object, label: str, path: str) -> str:
    return "|".join(str(part) for part in (task, backend, label, path))


def _artifact_kind(task: object, label: str) -> str:
    normalized = label.replace("-", "_")
    known = {
        "abundance": "taxonomy_table",
        "base_transition_plot": "base_transition_visualization",
        "base_transition_stats": "base_transition_stats",
        "bowtie2out": "alignment_index_output",
        "contig_map": "contig_map",
        "contigs": "contigs",
        "denoising_stats": "denoising_stats",
        "mags": "mag_bins",
        "output_dir": f"{task}_output_dir",
        "per_read": "per_read_taxonomy",
        "profile": "taxonomy_profile",
        "relative_abundance": "taxonomy_table",
        "representative_sequences": "representative_sequences",
        "stats": "method_stats",
        "table": "feature_table",
        "taxonomy": "taxonomy_table",
        "unbinned_contigs": "unbinned_contigs",
        "vector": "alpha_diversity",
        "distance_matrix": "beta_diversity",
    }
    return known.get(normalized, normalized)


def _format_from_path(path: str) -> str:
    suffixes = Path(path).suffixes
    if suffixes[-2:] == [".fastq", ".gz"]:
        return "fastq.gz"
    if suffixes[-2:] == [".fq", ".gz"]:
        return "fastq.gz"
    if not suffixes:
        return "directory"
    return suffixes[-1].lstrip(".")


def _utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
