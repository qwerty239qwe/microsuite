from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd

from microsuite.diversity.alpha import alpha_diversity as _alpha_diversity
from microsuite.diversity.beta import beta_diversity as _beta_diversity
from microsuite.diversity.ecology import (
    beta_significance as _beta_significance,
)
from microsuite.diversity.ecology import (
    beta_turnover as _beta_turnover,
)
from microsuite.diversity.ecology import (
    constrained_ordination as _constrained_ordination,
)
from microsuite.diversity.ecology import (
    gamma_diversity as _gamma_diversity,
)
from microsuite.diversity.ecology import (
    mantel_test as _mantel_test,
)
from microsuite.diversity.ecology import (
    taxa_turnover as _taxa_turnover,
)
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


def alpha_diversity(
    adata: ad.AnnData, *, metric: str = "shannon", tree: str | Path | None = None
) -> pd.DataFrame:
    return _alpha_diversity(adata, metric, tree=tree)


def beta_diversity(adata: ad.AnnData, *, metric: str = "bray-curtis") -> pd.DataFrame:
    return _beta_diversity(adata, metric)


def pcoa(distance_matrix: pd.DataFrame, *, dimensions: int = 3) -> pd.DataFrame:
    return _pcoa(distance_matrix, dimensions=dimensions)


def beta_significance(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    column: str,
    method: str = "permanova",
    permutations: int = 999,
    seed: int = 0,
) -> pd.DataFrame:
    return _beta_significance(
        distance_matrix,
        metadata,
        column=column,
        method=method,
        permutations=permutations,
        seed=seed,
    )


def mantel_test(
    matrix_a: pd.DataFrame,
    matrix_b: pd.DataFrame,
    *,
    method: str = "pearson",
    permutations: int = 999,
    seed: int = 0,
) -> pd.DataFrame:
    return _mantel_test(
        matrix_a,
        matrix_b,
        method=method,
        permutations=permutations,
        seed=seed,
    )


def gamma_diversity(
    adata: ad.AnnData, *, group: str, metric: str = "observed_features"
) -> pd.DataFrame:
    return _gamma_diversity(adata, group=group, metric=metric)


def beta_turnover(adata: ad.AnnData, *, level: str | None = None) -> pd.DataFrame:
    return _beta_turnover(adata, level=level)


def taxa_turnover(adata: ad.AnnData, *, group: str, level: str | None = None) -> pd.DataFrame:
    return _taxa_turnover(adata, group=group, level=level)


def constrained_ordination(
    adata: ad.AnnData,
    *,
    constraints: list[str],
    method: str = "rda",
    dimensions: int = 2,
) -> pd.DataFrame:
    return _constrained_ordination(
        adata,
        constraints=constraints,
        method=method,
        dimensions=dimensions,
    )
