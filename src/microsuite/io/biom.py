from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.metadata import read_indexed_tsv
from microsuite.io.taxonomy import add_taxonomy_levels, normalize_taxonomy_columns
from microsuite.io.tsv import read_tsv


def read_biom(
    table_path: Path, metadata_path: Path, taxonomy_path: Path | None = None
) -> ad.AnnData:
    try:
        biom_module = importlib.import_module("biom")
    except ImportError as exc:
        raise MicrobiomeSuiteError(
            "BIOM import requires the optional dependency 'biom-format'. "
            "Install with: uv sync --extra biom"
        ) from exc
    load_table = cast(Any, biom_module).load_table

    biom_table = load_table(str(table_path))
    sample_ids = [str(value) for value in biom_table.ids(axis="sample")]
    feature_ids = [str(value) for value in biom_table.ids(axis="observation")]
    matrix = biom_table.matrix_data.toarray().T.astype(np.float64)

    metadata = read_indexed_tsv(metadata_path, index_name="sample")
    missing = [sample for sample in sample_ids if sample not in metadata.index]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata is missing samples from BIOM table: {missing[:5]}")
    obs = metadata.loc[sample_ids].copy()
    var = pd.DataFrame(index=pd.Index(feature_ids, name="feature_id"))

    if taxonomy_path is not None:
        taxonomy = read_indexed_tsv(taxonomy_path, index_name="feature")
        taxonomy = taxonomy.rename(columns=normalize_taxonomy_columns(taxonomy.columns))
        var = var.join(taxonomy.reindex(var.index))
    else:
        taxonomy_values = []
        for feature_id in feature_ids:
            feature_metadata = biom_table.metadata(feature_id, axis="observation") or {}
            taxonomy = feature_metadata.get("taxonomy", "")
            if isinstance(taxonomy, list):
                taxonomy = "; ".join(str(part) for part in taxonomy)
            taxonomy_values.append(str(taxonomy))
        if any(taxonomy_values):
            var["taxonomy"] = taxonomy_values

    var = add_taxonomy_levels(var)
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    adata.uns["microsuite"] = {
        "version": "0.1.0",
        "importer": "biom",
        "table": str(table_path),
        "metadata": str(metadata_path),
        "taxonomy": str(taxonomy_path) if taxonomy_path else None,
    }
    return adata


def read_biom_or_tsv(
    table_path: Path, metadata_path: Path, taxonomy_path: Path | None = None
) -> ad.AnnData:
    if table_path.suffix.lower() in {".tsv", ".txt"}:
        return read_tsv(table_path, metadata_path, taxonomy_path)
    return read_biom(table_path, metadata_path, taxonomy_path)
