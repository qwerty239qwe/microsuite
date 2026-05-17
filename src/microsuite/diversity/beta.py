from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts

BETA_METRICS = {"bray-curtis", "jaccard"}


def beta_diversity(adata: ad.AnnData, metric: str) -> pd.DataFrame:
    metric = metric.lower().replace("_", "-")
    if metric not in BETA_METRICS:
        raise MicrobiomeSuiteError(
            f"Unsupported beta metric '{metric}'. Choose one of: {sorted(BETA_METRICS)}"
        )
    counts = dense_counts(adata)
    if metric == "bray-curtis":
        distances = _bray_curtis(counts)
    else:
        distances = _jaccard(counts)
    sample_ids = adata.obs_names.astype(str)
    return pd.DataFrame(distances, index=sample_ids, columns=sample_ids)


def _bray_curtis(counts: np.ndarray) -> np.ndarray:
    n_samples = counts.shape[0]
    distances = np.zeros((n_samples, n_samples), dtype=np.float64)
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            numerator = np.abs(counts[i] - counts[j]).sum()
            denominator = (counts[i] + counts[j]).sum()
            distance = 0.0 if denominator == 0 else numerator / denominator
            distances[i, j] = distances[j, i] = distance
    return distances


def _jaccard(counts: np.ndarray) -> np.ndarray:
    present = counts > 0
    n_samples = present.shape[0]
    distances = np.zeros((n_samples, n_samples), dtype=np.float64)
    for i in range(n_samples):
        for j in range(i + 1, n_samples):
            union = np.logical_or(present[i], present[j]).sum()
            intersection = np.logical_and(present[i], present[j]).sum()
            distance = 0.0 if union == 0 else 1.0 - (intersection / union)
            distances[i, j] = distances[j, i] = distance
    return distances
