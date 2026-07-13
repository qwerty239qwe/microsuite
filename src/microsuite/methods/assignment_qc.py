from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts
from microsuite.io.taxonomy import LEVELS

POOLED_LABEL = "_all_samples_"
SUMMARY_COLUMNS = [
    "sample",
    "rank",
    "assigned_features",
    "unassigned_features",
    "assigned_reads",
    "unassigned_reads",
    "assigned_feature_frac",
    "assigned_read_frac",
]


def _ranks(adata: ad.AnnData) -> list[str]:
    ranks = [r for r in LEVELS if r in adata.var.columns]
    if not ranks:
        raise MicrobiomeSuiteError(
            "Table has no taxonomy rank columns; run tax_classify or import with --taxonomy first."
        )
    return ranks


def _assigned_mask(series: pd.Series) -> np.ndarray:
    """Return a boolean mask of "assigned" entries.

    Both real NaN/NA values and empty or whitespace-only strings are treated
    as unassigned; anything else (a real taxon label) is assigned.
    """
    values = series.astype("object")
    empty = values.isna() | (values.fillna("").astype(str).str.strip() == "")
    return (~empty).to_numpy()


def _assigned_masks(adata: ad.AnnData, ranks: list[str]) -> dict[str, np.ndarray]:
    return {r: _assigned_mask(pd.Series(adata.var[r])) for r in ranks}


def _row(
    sample: str,
    present: np.ndarray,
    reads: np.ndarray,
    assigned: np.ndarray,
    rank: str,
) -> list:
    af = int((present & assigned).sum())
    uf = int((present & ~assigned).sum())
    ar = float(reads[assigned].sum())
    ur = float(reads[~assigned].sum())
    ff = af / (af + uf) if (af + uf) else 0.0
    rf = ar / (ar + ur) if (ar + ur) else 0.0
    return [sample, rank, af, uf, ar, ur, round(ff, 6), round(rf, 6)]


def summarize_assignment(adata: ad.AnnData) -> pd.DataFrame:
    ranks = _ranks(adata)
    counts = dense_counts(adata)  # samples x features
    assigned = _assigned_masks(adata, ranks)
    rows: list[list] = []
    for i, sample in enumerate([str(s) for s in adata.obs_names]):
        reads = counts[i]
        present = reads > 0
        for rank in ranks:
            rows.append(_row(sample, present, reads, assigned[rank], rank))
    pooled_reads = counts.sum(axis=0)
    pooled_present = pooled_reads > 0
    for rank in ranks:
        rows.append(_row(POOLED_LABEL, pooled_present, pooled_reads, assigned[rank], rank))
    return pd.DataFrame(rows, columns=pd.Index(SUMMARY_COLUMNS))


def write_assignment_summary(summary: pd.DataFrame, out_path: Path) -> Path:
    summary.to_csv(out_path, sep="\t", index=False)
    return out_path


def deepest_rank_distribution(adata: ad.AnnData) -> pd.Series:
    ranks = _ranks(adata)
    deepest = pd.Series("Unassigned", index=adata.var.index)
    for rank in ranks:  # shallow -> deep; later ranks overwrite
        mask = _assigned_mask(pd.Series(adata.var[rank]))
        deepest[mask] = rank
    return deepest.value_counts().reindex([*ranks, "Unassigned"], fill_value=0)
