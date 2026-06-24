from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import anndata as ad

from microsuite._errors import MicrobiomeSuiteError
from microsuite.api.table import read_table, write_table

AssetKind = Literal[
    "fasta",
    "fastq",
    "paired_fastq",
    "qza",
    "tree",
    "classifier",
    "reference_db",
]

SUPPORTED_ASSET_FORMATS: dict[AssetKind, tuple[str, ...]] = {
    "fasta": (".fa", ".fasta", ".fna", ".fas"),
    "fastq": (".fq", ".fastq", ".fq.gz", ".fastq.gz"),
    "paired_fastq": (".fq", ".fastq", ".fq.gz", ".fastq.gz"),
    "qza": (".qza",),
    "tree": (".nwk", ".newick", ".tree", ".qza"),
    "classifier": (".qza", ".pkl", ".pickle", ".joblib"),
    "reference_db": (),
}


@dataclass
class DataAsset:
    """Path-backed dataset input or result that should not be forced into AnnData."""

    name: str
    kind: AssetKind
    path: Path
    format: str
    sample_id: str | None = None
    paired_path: Path | None = None


class MicrobiomeDataset:
    """Stateful SDK object for table analysis plus registered raw data assets."""

    def __init__(
        self,
        adata: ad.AnnData | None = None,
        *,
        assets: dict[str, DataAsset] | None = None,
    ) -> None:
        self.adata = adata
        self.assets = dict(assets or {})

    @classmethod
    def read(cls, path: str | Path) -> MicrobiomeDataset:
        return cls(read_table(path))

    def write(self, path: str | Path) -> None:
        if self.adata is None:
            raise MicrobiomeSuiteError("MicrobiomeDataset has no feature table to write.")
        write_table(self.adata, path)

    def add_asset(
        self,
        name: str,
        path: str | Path,
        *,
        kind: AssetKind,
        sample_id: str | None = None,
        paired_path: str | Path | None = None,
        require_exists: bool = True,
    ) -> DataAsset:
        if name in self.assets:
            raise MicrobiomeSuiteError(f"Dataset asset already exists: {name}")

        asset_path = Path(path)
        mate_path = Path(paired_path) if paired_path is not None else None
        if require_exists:
            _require_existing_path(asset_path)
            if mate_path is not None:
                _require_existing_path(mate_path)

        asset_format = _asset_format(asset_path)
        _validate_asset_format(kind, asset_format)
        if kind == "paired_fastq":
            if mate_path is None:
                raise MicrobiomeSuiteError("paired_fastq assets require paired_path.")
            mate_format = _asset_format(mate_path)
            _validate_asset_format(kind, mate_format)
            if mate_format != asset_format:
                raise MicrobiomeSuiteError(
                    "paired_fastq assets require r1 and r2 to use the same format."
                )

        asset = DataAsset(
            name=name,
            kind=kind,
            path=asset_path,
            paired_path=mate_path,
            format=asset_format,
            sample_id=sample_id,
        )
        self.assets[name] = asset
        return asset

    def add_reads(
        self,
        sample_id: str,
        *,
        r1: str | Path,
        r2: str | Path | None = None,
        name: str | None = None,
        require_exists: bool = True,
    ) -> DataAsset:
        asset_name = name or sample_id
        if r2 is None:
            return self.add_asset(
                asset_name,
                r1,
                kind="fastq",
                sample_id=sample_id,
                require_exists=require_exists,
            )
        return self.add_asset(
            asset_name,
            r1,
            kind="paired_fastq",
            sample_id=sample_id,
            paired_path=r2,
            require_exists=require_exists,
        )

    def add_sequences(
        self,
        name: str,
        path: str | Path,
        *,
        require_exists: bool = True,
    ) -> DataAsset:
        return self.add_asset(name, path, kind="fasta", require_exists=require_exists)

    def add_tree(
        self,
        name: str,
        path: str | Path,
        *,
        require_exists: bool = True,
    ) -> DataAsset:
        return self.add_asset(name, path, kind="tree", require_exists=require_exists)

    def add_artifact(
        self,
        name: str,
        path: str | Path,
        *,
        require_exists: bool = True,
    ) -> DataAsset:
        return self.add_asset(name, path, kind="qza", require_exists=require_exists)

    def add_reference(
        self,
        name: str,
        path: str | Path,
        *,
        kind: AssetKind = "reference_db",
        require_exists: bool = True,
    ) -> DataAsset:
        return self.add_asset(name, path, kind=kind, require_exists=require_exists)

    def list_assets(self, *, kind: AssetKind | None = None) -> list[DataAsset]:
        assets = list(self.assets.values())
        if kind is None:
            return assets
        return [asset for asset in assets if asset.kind == kind]

    def get_asset(self, name: str) -> DataAsset:
        try:
            return self.assets[name]
        except KeyError as exc:
            raise MicrobiomeSuiteError(f"Unknown dataset asset: {name}") from exc

    def validate_assets(self) -> None:
        for asset in self.assets.values():
            _require_existing_path(asset.path)
            _validate_asset_format(asset.kind, asset.format)
            if asset.paired_path is not None:
                _require_existing_path(asset.paired_path)
                _validate_asset_format(asset.kind, _asset_format(asset.paired_path))


def _asset_format(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".fastq.gz"):
        return ".fastq.gz"
    if name.endswith(".fq.gz"):
        return ".fq.gz"
    return path.suffix.lower()


def _validate_asset_format(kind: AssetKind, asset_format: str) -> None:
    if kind == "reference_db":
        return
    supported = SUPPORTED_ASSET_FORMATS[kind]
    if asset_format not in supported:
        expected = ", ".join(supported)
        raise MicrobiomeSuiteError(
            f"{kind} asset expects one of {expected}; got {asset_format or 'no suffix'}."
        )


def _require_existing_path(path: Path) -> None:
    if not path.exists():
        raise MicrobiomeSuiteError(f"Dataset asset does not exist: {path}")
