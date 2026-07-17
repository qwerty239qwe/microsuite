from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.metadata import (
    sha256_file,
    validate_reads_manifest,
    validate_run_bundle,
    validate_run_manifest,
    validate_workflow,
    write_reads_manifest,
    write_resolved_config,
    write_run_manifest,
    write_workflow,
)
from microsuite.metadata.schemas import (
    READS_MANIFEST_VERSION,
    RUN_MANIFEST_VERSION,
    WORKFLOW_VERSION,
    published_schema_path,
)


def _workflow() -> dict:
    return {
        "schema_version": WORKFLOW_VERSION,
        "workflow_id": "oral-amplicon",
        "workflow_version": "1.0.0",
        "workflow_hash": "sha256:abc",
        "generated_at": "2026-07-17T00:00:00Z",
        "producer": {"name": "microsuite", "version": "0.2.0.dev0"},
        "engine": {"name": "nextflow", "version": "25.10.4"},
        "nodes": [
            {
                "id": "cutadapt",
                "stage_type": "primer_trimming",
                "process_name": "ORAL:READS:CUTADAPT",
                "method": "cutadapt",
                "branch": "reads",
                "rank": None,
                "inputs": ["raw_reads"],
                "outputs": ["trimmed_reads"],
            }
        ],
        "edges": [],
        "branches": [{"id": "reads", "enabled": True, "node_ids": ["cutadapt"]}],
        "methods": [{"id": "cutadapt", "name": "cutadapt", "version": "4.9"}],
        "intermediate_files": [
            {
                "id": "trimmed_reads",
                "label": "Primer-trimmed reads",
                "kind": "reads",
                "format": "fastq.gz",
                "path": "reads/trimmed",
                "producer_node": "cutadapt",
                "consumer_nodes": [],
                "retained": True,
            }
        ],
        "native_dag": {"path": "provenance/nextflow/dag.dot", "sha256": None},
    }


def _reads() -> dict:
    return {
        "schema_version": READS_MANIFEST_VERSION,
        "dataset_id": "ERP120510",
        "generated_at": "2026-07-17T00:00:00Z",
        "producer": {"name": "microsuite", "version": "0.2.0.dev0"},
        "layout": "PE",
        "samples": [
            {
                "sample_id": "S1",
                "layout": "PE",
                "files": [
                    {"path": "/data/S1_R1.fastq.gz", "read": "R1", "external": True},
                    {"path": "/data/S1_R2.fastq.gz", "read": "R2", "external": True},
                ],
                "metadata": {"site_code": "SubG"},
            }
        ],
        "selection": {"column": "site_code", "value": "SubG"},
        "source": {"input_dir": "/data"},
    }


def _manifest(run_dir: Path) -> dict:
    return {
        "schema_version": RUN_MANIFEST_VERSION,
        "run_id": "run-1",
        "dataset_id": "ERP120510",
        "analysis_id": "v4v5-subg",
        "status": "running",
        "created_at": "2026-07-17T00:00:00Z",
        "updated_at": "2026-07-17T00:00:00Z",
        "producer": {"name": "microsuite", "version": "0.2.0.dev0"},
        "workflow": {
            "path": "workflow.json",
            "schema_version": WORKFLOW_VERSION,
            "sha256": sha256_file(run_dir / "workflow.json"),
        },
        "resolved_config": {
            "path": "resolved_config.json",
            "schema_version": "resolved-config.v1",
            "sha256": sha256_file(run_dir / "resolved_config.json"),
        },
        "reads_manifest": {
            "path": "reads-manifest.v1.json",
            "schema_version": READS_MANIFEST_VERSION,
            "sha256": sha256_file(run_dir / "reads-manifest.v1.json"),
        },
        "progress": {
            "total": 1,
            "completed": 0,
            "failed": 0,
            "running": 0,
            "pending": 1,
            "fraction": 0.0,
        },
        "stages": [
            {
                "node_id": "cutadapt",
                "status": "pending",
                "expected_tasks": 1,
                "completed_tasks": 0,
                "failed_tasks": 0,
                "result_files": [],
            }
        ],
        "artifacts": [],
        "summary_metrics": {},
        "provenance": {"project_revision": "abc123"},
        "error": None,
    }


def test_workflow_and_reads_cross_field_validation() -> None:
    workflow = _workflow()
    workflow["edges"] = [{"source": "cutadapt", "target": "missing"}]
    assert any("unknown node" in error for error in validate_workflow(workflow))

    reads = _reads()
    reads["samples"][0]["files"].pop()
    assert any("PE layout" in error for error in validate_reads_manifest(reads))


def test_run_manifest_progress_and_path_validation(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    write_workflow(run_dir / "workflow.json", _workflow())
    write_resolved_config(run_dir, {"project": {"accession": "ERP120510"}})
    write_reads_manifest(run_dir / "reads-manifest.v1.json", _reads())
    manifest = _manifest(run_dir)
    manifest["progress"]["fraction"] = 1.0
    assert any("progress.fraction" in error for error in validate_run_manifest(manifest))
    manifest = _manifest(run_dir)
    manifest["workflow"]["path"] = "../workflow.json"
    assert any("safe relative path" in error for error in validate_run_manifest(manifest))


def test_writers_and_bundle_validator(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    write_workflow(run_dir / "workflow.json", _workflow())
    write_resolved_config(run_dir, {"project": {"accession": "ERP120510"}})
    write_reads_manifest(run_dir / "reads-manifest.v1.json", _reads())
    write_run_manifest(run_dir / "run_manifest.json", _manifest(run_dir))

    assert validate_run_bundle(run_dir) == []
    loaded = json.loads((run_dir / "run_manifest.json").read_text(encoding="utf-8"))
    assert loaded["schema_version"] == RUN_MANIFEST_VERSION


def test_bundle_validator_detects_checksum_mismatch(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    write_workflow(run_dir / "workflow.json", _workflow())
    write_resolved_config(run_dir, {})
    write_reads_manifest(run_dir / "reads-manifest.v1.json", _reads())
    manifest = _manifest(run_dir)
    manifest["workflow"]["sha256"] = "0" * 64
    write_run_manifest(run_dir / "run_manifest.json", manifest)
    assert any("checksum mismatch" in error for error in validate_run_bundle(run_dir))


def test_writer_refuses_invalid_document(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="invalid workflow.v1"):
        write_workflow(tmp_path / "workflow.json", {"workflow_id": "missing-fields"})


@pytest.mark.parametrize(
    ("schema_name", "payload"),
    [
        (WORKFLOW_VERSION, _workflow()),
        (READS_MANIFEST_VERSION, _reads()),
    ],
)
def test_published_schemas_accept_writer_documents(
    tmp_path: Path, schema_name: str, payload: dict
) -> None:
    jsonschema = pytest.importorskip("jsonschema")
    target = tmp_path / f"{schema_name}.json"
    if schema_name == WORKFLOW_VERSION:
        write_workflow(target, payload)
    else:
        write_reads_manifest(target, payload)
    schema = json.loads(published_schema_path(schema_name).read_text(encoding="utf-8"))
    cls = jsonschema.validators.validator_for(schema)
    cls.check_schema(schema)
    cls(schema, format_checker=cls.FORMAT_CHECKER).validate(
        json.loads(target.read_text(encoding="utf-8"))
    )
