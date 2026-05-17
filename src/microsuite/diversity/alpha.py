from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts

ALPHA_METRICS = {"observed_features", "shannon", "simpson", "pielou"}


def alpha_diversity(adata: ad.AnnData, metric: str) -> pd.DataFrame:
    metric = metric.lower().replace("_", "-")
    canonical = metric.replace("-", "_")
    if canonical not in ALPHA_METRICS:
        raise MicrobiomeSuiteError(
            f"Unsupported alpha metric '{metric}'. Choose one of: {sorted(ALPHA_METRICS)}"
        )
    counts = dense_counts(adata)
    values = _compute(counts, canonical)
    return pd.DataFrame({"sample_id": adata.obs_names.astype(str), canonical: values})


def _compute(counts: np.ndarray, metric: str) -> np.ndarray:
    totals = counts.sum(axis=1)
    observed = (counts > 0).sum(axis=1).astype(np.float64)
    proportions = np.zeros_like(counts, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        np.divide(counts, totals[:, None], out=proportions, where=totals[:, None] > 0)

    if metric == "observed_features":
        return observed
    if metric == "shannon":
        with np.errstate(divide="ignore", invalid="ignore"):
            terms = np.where(proportions > 0, proportions * np.log(proportions), 0.0)
        return -terms.sum(axis=1)
    if metric == "simpson":
        return 1.0 - np.square(proportions).sum(axis=1)
    if metric == "pielou":
        shannon = _compute(counts, "shannon")
        with np.errstate(divide="ignore", invalid="ignore"):
            evenness = np.divide(shannon, np.log(observed), where=observed > 1)
        evenness[observed <= 1] = 0.0
        return evenness
    raise AssertionError(metric)
