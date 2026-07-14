from __future__ import annotations

from pathlib import Path

import pytest

from microsuite.metadata.models import Artifact, ArtifactCount, ProvenanceFile, StageError


def test_artifact_accepts_absolute_path() -> None:
    art = Artifact(label="table", path="/abs/table.tsv", format="tsv")
    assert art.required is True
    assert art.external is False
    assert art.count is None


def test_artifact_accepts_pathlib_absolute() -> None:
    art = Artifact(label="table", path=Path("/abs/table.tsv"))
    assert Path(art.path).is_absolute()


def test_artifact_rejects_relative_path() -> None:
    with pytest.raises(ValueError, match="must be absolute"):
        Artifact(label="table", path="table.tsv")


def test_provenance_rejects_relative_path() -> None:
    with pytest.raises(ValueError, match="must be absolute"):
        ProvenanceFile(kind="dada2_manifest", path="manifest.json")


def test_provenance_accepts_absolute() -> None:
    prov = ProvenanceFile(kind="dada2_manifest", path="/abs/manifest.json")
    assert prov.required is True


def test_artifact_count_and_stage_error() -> None:
    count = ArtifactCount(value=1842, unit="features")
    assert count.value == 1842 and count.unit == "features"
    err = StageError(type="ValueError", message="boom")
    assert err.type == "ValueError" and err.message == "boom"
