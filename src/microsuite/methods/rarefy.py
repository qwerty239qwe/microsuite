from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity._matrix import dense_counts
from microsuite.io.h5ad import read_h5ad, write_h5ad
from microsuite.methods._dispatch import require_backend

SUPPORTED_BACKENDS = ("native",)


def rarefy(
    *,
    backend: str,
    table: Path,
    output: Path,
    depth: int,
    seed: int = 0,
    force: bool = False,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "rarefy")
    result = rarefy_native(read_h5ad(ensure_input(table)), depth=depth, seed=seed)
    write_h5ad(result, prepare_output(output, force=force))


def rarefy_native(adata: ad.AnnData, *, depth: int, seed: int = 0) -> ad.AnnData:
    if depth <= 0:
        raise MicrobiomeSuiteError("Rarefaction depth must be greater than zero.")
    counts = dense_counts(adata)
    totals = counts.sum(axis=1)
    low_depth = adata.obs_names[totals < depth].astype(str).tolist()
    if low_depth:
        raise MicrobiomeSuiteError(
            f"Samples below rarefaction depth {depth}: {', '.join(low_depth[:5])}"
        )

    rng = np.random.default_rng(seed)
    rarefied = np.zeros_like(counts, dtype=np.float64)
    n_features = counts.shape[1]
    for index, row in enumerate(counts):
        # Subsample WITHOUT replacement: expand the row into one entry per read,
        # shuffle, keep the first `depth`, then re-bin. Matches QIIME 2 / scikit-bio
        # rarefaction (a with-replacement draw would allow a feature to exceed its
        # original count and would not return the table unchanged at depth == total).
        reads = np.repeat(np.arange(n_features), row.astype(np.int64))
        rng.shuffle(reads)
        rarefied[index] = np.bincount(reads[:depth], minlength=n_features)

    result = cast(Any, adata).copy()
    result.X = rarefied
    result.uns["microsuite_rarefy"] = {
        "backend": "native",
        "depth": depth,
        "seed": seed,
    }
    return result
