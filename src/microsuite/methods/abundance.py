from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity._matrix import dense_counts
from microsuite.io.h5ad import read_h5ad
from microsuite.methods._dispatch import require_backend

SUPPORTED_BACKENDS = ("native",)


def abundance(
    *,
    backend: str,
    table: Path,
    output: Path,
    level: str,
    relative: bool = True,
    force: bool = False,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "abundance")
    result = abundance_native(read_h5ad(ensure_input(table)), level=level, relative=relative)
    result.to_csv(prepare_output(output, force=force), sep="\t")


def abundance_native(adata: ad.AnnData, *, level: str, relative: bool = True) -> pd.DataFrame:
    var = pd.DataFrame(adata.var)
    if level not in var.columns:
        raise MicrobiomeSuiteError(f"Taxonomy level not found in table: {level}")
    taxa = var[level].astype(str).replace("", "Unclassified")
    frame = pd.DataFrame(dense_counts(adata), index=adata.obs_names.astype(str), columns=taxa)
    grouped = frame.T.groupby(level=0).sum().T
    if relative:
        totals = grouped.sum(axis=1).replace(0, np.nan)
        grouped = grouped.div(totals, axis=0).fillna(0.0)
    grouped.index.name = "sample_id"
    return grouped
