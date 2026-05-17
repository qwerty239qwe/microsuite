from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from zipfile import ZipFile

import anndata as ad

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.biom import read_biom_or_tsv
from microsuite.qiime2.artifact import inspect_artifact


def read_qza(
    artifact_path: Path,
    metadata_path: Path,
    taxonomy_artifact_path: Path | None = None,
) -> ad.AnnData:
    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        table_path = _extract_table_artifact(artifact_path, temp / "table")
        taxonomy_path = (
            _extract_taxonomy_artifact(taxonomy_artifact_path, temp / "taxonomy")
            if taxonomy_artifact_path
            else None
        )
        adata = read_biom_or_tsv(table_path, metadata_path, taxonomy_path)
        table_info = inspect_artifact(artifact_path)
        adata.uns["qiime2_table_artifact"] = {
            "uuid": table_info.uuid,
            "type": table_info.artifact_type,
            "format": table_info.format,
            "framework_version": table_info.framework_version,
            "archive_version": table_info.archive_version,
            "data_files": table_info.data_files,
        }
        if taxonomy_artifact_path:
            taxonomy_info = inspect_artifact(taxonomy_artifact_path)
            adata.uns["qiime2_taxonomy_artifact"] = {
                "uuid": taxonomy_info.uuid,
                "type": taxonomy_info.artifact_type,
                "format": taxonomy_info.format,
                "framework_version": taxonomy_info.framework_version,
                "archive_version": taxonomy_info.archive_version,
                "data_files": taxonomy_info.data_files,
            }
        return adata


def _extract_table_artifact(artifact_path: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    with ZipFile(artifact_path) as archive:
        names = archive.namelist()
        table_names = [
            name
            for name in names
            if name.endswith("data/feature-table.biom")
            or name.endswith("data/table.biom")
            or name.endswith("data/feature-frequency.biom")
            or name.endswith("data/feature-table.tsv")
            or name.endswith("data/table.tsv")
            or name.endswith("data/data.tsv")
        ]
        if not table_names:
            raise MicrobiomeSuiteError(
                "QZA artifact does not contain a supported feature table "
                "(expected data/feature-table.biom or data/feature-table.tsv)."
            )
        member = table_names[0]
        suffix = Path(member).suffix
        target = output_dir / f"feature-table{suffix}"
        target.write_bytes(archive.read(member))
        return target


def _extract_taxonomy_artifact(artifact_path: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    with ZipFile(artifact_path) as archive:
        names = archive.namelist()
        taxonomy_names = [
            name
            for name in names
            if name.endswith("data/taxonomy.tsv")
            or name.endswith("data/feature-taxonomy.tsv")
            or name.endswith("data/data.tsv")
        ]
        if not taxonomy_names:
            raise MicrobiomeSuiteError("Taxonomy QZA artifact does not contain data/taxonomy.tsv.")
        target = output_dir / "taxonomy.tsv"
        target.write_bytes(archive.read(taxonomy_names[0]))
        return target
