from __future__ import annotations

from datetime import UTC, datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from microsuite import __version__
from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.metadata import read_indexed_tsv
from microsuite.io.taxonomy import add_taxonomy_levels, normalize_taxonomy_columns


def read_tsv(
    table_path: Path, metadata_path: Path, taxonomy_path: Path | None = None
) -> ad.AnnData:
    table = pd.read_csv(table_path, sep="\t")
    if table.empty or table.shape[1] < 2:
        raise MicrobiomeSuiteError(
            "Feature table must have a feature ID column and sample columns."
        )

    feature_col = table.columns[0]
    table = table.set_index(feature_col)
    table.index = table.index.astype(str)
    if table.index.has_duplicates:
        raise MicrobiomeSuiteError("Feature table contains duplicate feature IDs.")

    counts = table.apply(pd.to_numeric, errors="raise")
    counts.columns = counts.columns.astype(str)

    metadata = read_indexed_tsv(metadata_path, index_name="sample")
    missing = [sample for sample in counts.columns if sample not in metadata.index]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata is missing samples from table: {missing[:5]}")
    metadata = metadata.loc[counts.columns].copy()

    var = pd.DataFrame(index=counts.index.astype(str))
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
