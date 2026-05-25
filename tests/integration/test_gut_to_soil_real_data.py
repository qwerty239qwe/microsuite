from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
import pytest

from microsuite.data.gut_to_soil import fetch_gut_to_soil
from microsuite.diversity.alpha import alpha_diversity
from microsuite.diversity.beta import beta_diversity
from microsuite.io.h5ad import write_h5ad
from microsuite.io.qza import read_qza
from microsuite.methods.abundance import abundance_native
from microsuite.methods.ml_longitudinal import ml_classify_native
from microsuite.methods.network import correlation_network
from microsuite.methods.normalize import normalize_native
from microsuite.ordination.pcoa import pcoa

pytestmark = pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_REAL_DATA_TESTS") != "1",
    reason="set MICROSUITE_RUN_REAL_DATA_TESTS=1 to download and run real-data tests",
)


def test_gut_to_soil_real_data_native_methods(tmp_path: Path) -> None:
    pytest.importorskip("biom")
    data_dir = tmp_path / "gut-to-soil"
    fetch_gut_to_soil(data_dir)

    adata = read_qza(
        data_dir / "asv-table-ms2.qza",
        data_dir / "sample-metadata.tsv",
        data_dir / "taxonomy.qza",
    )
    table_path = tmp_path / "gut-to-soil.h5ad"
    write_h5ad(adata, table_path)

    assert adata.n_obs > 10
    assert adata.n_vars > 10
    assert "qiime2_table_artifact" in adata.uns

    normalized = normalize_native(adata, method="relative")
    assert normalized.shape == adata.shape

    alpha = alpha_diversity(adata, "shannon")
    beta = beta_diversity(adata, "bray-curtis")
    coords = pcoa(beta, dimensions=2)

    assert alpha.shape[0] == adata.n_obs
    assert beta.shape == (adata.n_obs, adata.n_obs)
    assert coords.shape[0] == adata.n_obs

    genus = abundance_native(adata, level="genus", relative=True)
    edges = correlation_network(
        adata,
        method="spearman",
        transform="relative",
        min_abs_weight=0.6,
        min_prevalence=0.05,
        top_n=25,
    )

    assert genus.shape[0] == adata.n_obs
    assert {"source", "target", "weight"}.issubset(edges.columns)

    target, keep_values = _categorical_metadata_column(pd.DataFrame(adata.obs))
    ml_adata = adata[adata.obs[target].astype(str).isin(keep_values)].copy()
    predictions, importance = ml_classify_native(
        ml_adata,
        backend="randomforest",
        target=target,
        test_fraction=0.25,
        n_estimators=25,
        seed=7,
    )

    assert not predictions.empty
    assert {"truth", "prediction", "correct"}.issubset(predictions.columns)
    assert not importance.empty


def _categorical_metadata_column(metadata: pd.DataFrame) -> tuple[str, list[str]]:
    for column in ["sample-type", "SampleType", "body-site", "BodySite"]:
        if column in metadata.columns and metadata[column].astype(str).nunique() >= 2:
            counts = metadata[column].astype(str).value_counts()
            if (counts >= 2).sum() >= 2:
                return column, counts[counts >= 2].index.astype(str).tolist()[:3]
    for column in metadata.columns:
        values = metadata[column].astype(str)
        counts = values.value_counts()
        if 2 <= values.nunique() <= max(10, len(values) // 2) and (counts >= 2).sum() >= 2:
            return column, counts[counts >= 2].index.astype(str).tolist()[:3]
    raise AssertionError("No suitable categorical metadata column found for classification.")
