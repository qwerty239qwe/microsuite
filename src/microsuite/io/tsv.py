from __future__ import annotations

import warnings
from datetime import UTC, datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from microsuite import __version__
from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.metadata import read_indexed_tsv
from microsuite.io.taxonomy import LEVELS, add_taxonomy_levels, normalize_taxonomy_columns

_RESERVED_FEATURE_NAMES = set(LEVELS) | {"taxonomy", "taxon"}


def read_count_matrix(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    if table.empty or table.shape[1] < 2:
        raise MicrobiomeSuiteError(
            "Feature table must have a feature ID column and sample columns."
        )
    original = str(table.columns[0])
    counts = table.set_index(table.columns[0])
    counts.index = counts.index.astype(str)
    if counts.index.has_duplicates:
        raise MicrobiomeSuiteError("Feature table contains duplicate feature IDs.")
    counts = counts.apply(pd.to_numeric, errors="raise")
    counts.columns = counts.columns.astype(str)
    counts.index.name = "feature_id"
    if original.strip().lower() in _RESERVED_FEATURE_NAMES:
        warnings.warn(
            f"Renamed feature-ID column '{original}' to 'feature_id' to avoid a "
            "taxonomy-rank naming conflict; feature IDs are unchanged.",
            stacklevel=2,
        )
    return counts


def read_matrix_tsv(path: Path) -> ad.AnnData:
    counts = read_count_matrix(path)
    obs = pd.DataFrame(index=counts.columns)
    var = add_taxonomy_levels(pd.DataFrame(index=counts.index))
    adata = ad.AnnData(X=counts.T.to_numpy(dtype=np.float64), obs=obs, var=var)
    adata.uns["microsuite"] = {
        "version": __version__,
        "importer": "matrix-tsv",
        "table": str(path),
        "created_at": datetime.now(UTC).isoformat(),
    }
    return adata


def read_tsv(
    table_path: Path, metadata_path: Path, taxonomy_path: Path | None = None
) -> ad.AnnData:
    counts = read_count_matrix(table_path)

    metadata = read_indexed_tsv(metadata_path, index_name="sample")
    missing = [sample for sample in counts.columns if sample not in metadata.index]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata is missing samples from table: {missing[:5]}")
    metadata = metadata.loc[counts.columns].copy()

    var = pd.DataFrame(index=counts.index)
    if taxonomy_path is not None:
        taxonomy = read_indexed_tsv(taxonomy_path, index_name="feature")
        taxonomy = taxonomy.rename(columns=normalize_taxonomy_columns(taxonomy.columns))
        taxonomy = taxonomy.reindex(var.index)
        var = var.join(taxonomy)
    var = add_taxonomy_levels(var)

    adata = ad.AnnData(
        X=counts.T.to_numpy(dtype=np.float64),
        obs=metadata,
        var=var,
    )
    adata.uns["microsuite"] = {
        "version": __version__,
        "importer": "tsv",
        "table": str(table_path),
        "metadata": str(metadata_path),
        "taxonomy": str(taxonomy_path) if taxonomy_path else None,
        "created_at": datetime.now(UTC).isoformat(),
    }
    return adata
