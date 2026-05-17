from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import anndata as ad
import numpy as np

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity._matrix import dense_counts
from microsuite.io.h5ad import read_h5ad, write_h5ad

SUPPORTED_BACKENDS = ("native",)
NORMALIZE_METHODS = ("relative", "total-sum", "clr", "prevalence-filter")


def normalize(
    *,
    backend: str,
    method: str,
    table: Path,
    output: Path,
    target_sum: float = 1_000_000.0,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.1,
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend != "native":
        raise MicrobiomeSuiteError(
            f"Unsupported normalize backend '{backend}'. "
            f"Choose one of: {', '.join(SUPPORTED_BACKENDS)}"
        )
    adata = normalize_native(
        read_h5ad(ensure_input(table)),
        method=method,
        target_sum=target_sum,
        pseudocount=pseudocount,
        min_prevalence=min_prevalence,
    )
    write_h5ad(adata, prepare_output(output, force=force))


def normalize_native(
    adata: ad.AnnData,
    *,
    method: str,
    target_sum: float = 1_000_000.0,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.1,
) -> ad.AnnData:
    method = method.lower()
    if method not in NORMALIZE_METHODS:
        raise MicrobiomeSuiteError(
            f"Unsupported normalize method '{method}'. "
            f"Choose one of: {', '.join(NORMALIZE_METHODS)}"
        )
    counts = dense_counts(adata)
    result = cast(Any, adata).copy()

    if method == "relative":
        totals = counts.sum(axis=1)
        normalized = np.zeros_like(counts, dtype=np.float64)
        np.divide(counts, totals[:, None], out=normalized, where=totals[:, None] > 0)
        result.X = normalized
    elif method == "total-sum":
        totals = counts.sum(axis=1)
        normalized = np.zeros_like(counts, dtype=np.float64)
        np.divide(
            counts * target_sum,
            totals[:, None],
            out=normalized,
            where=totals[:, None] > 0,
        )
        result.X = normalized
    elif method == "clr":
        if pseudocount <= 0:
            raise MicrobiomeSuiteError("CLR pseudocount must be greater than zero.")
        logged = np.log(counts + pseudocount)
        result.X = logged - logged.mean(axis=1, keepdims=True)
    else:
        if not 0 <= min_prevalence <= 1:
            raise MicrobiomeSuiteError("min_prevalence must be between 0 and 1.")
        prevalence = (counts > 0).mean(axis=0)
        keep = prevalence >= min_prevalence
        result = result[:, keep].copy()

    result.uns["microsuite_normalize"] = {
        "backend": "native",
        "method": method,
        "target_sum": target_sum,
        "pseudocount": pseudocount,
        "min_prevalence": min_prevalence,
    }
    return result
