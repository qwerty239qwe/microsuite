from __future__ import annotations

import anndata as ad
import pandas as pd

from microsuite.diversity.alpha import alpha_diversity as _alpha_diversity
from microsuite.diversity.beta import beta_diversity as _beta_diversity
from microsuite.methods.abundance import abundance_native
from microsuite.methods.normalize import normalize_native
from microsuite.methods.rarefy import rarefy_native
from microsuite.methods.shared_taxa import shared_taxa_native
from microsuite.ordination.pcoa import pcoa as _pcoa


def normalize_table(
    adata: ad.AnnData,
    *,
    method: str = "relative",
    target_sum: float = 1_000_000.0,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.1,
) -> ad.AnnData:
    return normalize_native(
        adata,
        method=method,
        target_sum=target_sum,
        pseudocount=pseudocount,
        min_prevalence=min_prevalence,
    )


def abundance_table(
    adata: ad.AnnData, *, level: str = "genus", relative: bool = True
) -> pd.DataFrame:
    return abundance_native(adata, level=level, relative=relative)


def shared_taxa_table(adata: ad.AnnData, *, level: str, group: str) -> pd.DataFrame:
    return shared_taxa_native(adata, level=level, group=group)


def rarefy_table(adata: ad.AnnData, *, depth: int, seed: int = 0) -> ad.AnnData:
    return rarefy_native(adata, depth=depth, seed=seed)


def alpha_diversity(adata: ad.AnnData, *, metric: str = "shannon") -> pd.DataFrame:
    return _alpha_diversity(adata, metric)


def beta_diversity(adata: ad.AnnData, *, metric: str = "bray-curtis") -> pd.DataFrame:
    return _beta_diversity(adata, metric)


def pcoa(distance_matrix: pd.DataFrame, *, dimensions: int = 3) -> pd.DataFrame:
    return _pcoa(distance_matrix, dimensions=dimensions)
