from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile

from microsuite.io.qza import read_qza
from microsuite.qiime2.artifact import extract_data_payload, inspect_artifact

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def write_qiime_artifact(path: Path, *, payload: Path, payload_name: str) -> None:
    with ZipFile(path, "w") as archive:
        archive.writestr("uuid/VERSION", "QIIME 2\narchive: 5\nframework: 2024.10.0\n")
        archive.writestr(
            "uuid/metadata.yaml",
            "uuid: uuid\ntype: FeatureTable[Frequency]\nformat: BIOMV210DirFmt\n",
        )
        archive.write(payload, f"uuid/data/{payload_name}")


def test_qza_reads_minimal_tsv_artifact(tmp_path: Path) -> None:
    artifact = tmp_path / "table.qza"
    write_qiime_artifact(artifact, payload=FIXTURE / "table.tsv", payload_name="feature-table.tsv")

    adata = read_qza(artifact, FIXTURE / "metadata.tsv")

    assert adata.shape == (4, 4)
    assert list(adata.obs_names) == ["L1S8", "L1S57", "L1S76", "L2S155"]
    assert adata.uns["qiime2_table_artifact"]["uuid"] == "uuid"
    assert adata.uns["qiime2_table_artifact"]["type"] == "FeatureTable[Frequency]"


def test_qiime_artifact_inspect_and_extract(tmp_path: Path) -> None:
    artifact = tmp_path / "table.qza"
    write_qiime_artifact(artifact, payload=FIXTURE / "table.tsv", payload_name="data.tsv")

    info = inspect_artifact(artifact)
    extracted = extract_data_payload(artifact, tmp_path / "payload")

    assert info.uuid == "uuid"
    assert info.artifact_type == "FeatureTable[Frequency]"
    assert info.framework_version == "2024.10.0"
    assert info.archive_version == "5"
    assert info.data_files == ["data.tsv"]
    assert extracted == [tmp_path / "payload" / "data.tsv"]
    assert extracted[0].read_text(encoding="utf-8").startswith("feature_id")
