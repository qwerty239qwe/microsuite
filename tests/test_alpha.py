from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity.alpha import (
    alpha_diversity,
    available_alpha_metrics,
    chao1,
    dominance,
    goods_coverage,
    inv_simpson,
    margalef,
    menhinick,
    osd,
    pielou_e,
    simpson,
)
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def test_alpha_observed_features() -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    result = alpha_diversity(adata, "observed_features")

    assert result["observed_features"].tolist() == [2.0, 3.0, 3.0, 3.0]


def test_alpha_shannon_known_value() -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    result = alpha_diversity(adata, "shannon")

    expected = -((10 / 12) * np.log(10 / 12) + (2 / 12) * np.log(2 / 12))
    assert np.isclose(result.loc[0, "shannon"], expected)


def test_alpha_metric_registry_covers_scikit_bio_alpha_names() -> None:
    expected = {
        "ace",
        "chao1",
        "chao1_ci",
        "doubles",
        "faith_pd",
        "margalef",
        "menhinick",
        "michaelis_menten_fit",
        "observed_features",
        "osd",
        "singles",
        "sobs",
        "brillouin_d",
        "enspie",
        "fisher_alpha",
        "hill",
        "inv_simpson",
        "kempton_taylor_q",
        "phydiv",
        "renyi",
        "shannon",
        "simpson",
        "tsallis",
        "heip_e",
        "mcintosh_e",
        "pielou_e",
        "simpson_e",
        "berger_parker_d",
        "dominance",
        "gini_index",
        "mcintosh_d",
        "simpson_d",
        "strong",
        "esty_ci",
        "goods_coverage",
        "lladser_ci",
        "lladser_pe",
        "robbins",
    }

    assert expected.issubset(set(available_alpha_metrics()))


def test_alpha_scalar_formulas_on_known_counts() -> None:
    counts = np.array([4, 3, 2, 1, 0])

    assert osd(counts) == (4.0, 1.0, 1.0)
    assert np.isclose(chao1(counts), 4.0)
    assert np.isclose(dominance(counts), 0.3)
    assert np.isclose(simpson(counts), 0.7)
    assert np.isclose(inv_simpson(counts), 1 / 0.3)
    assert np.isclose(goods_coverage(counts), 0.9)
    assert np.isclose(margalef(counts), (4 - 1) / np.log(10))
    assert np.isclose(menhinick(counts), 4 / np.sqrt(10))
    proportions = counts[counts > 0] / 10
    assert np.isclose(pielou_e(counts), -np.sum(proportions * np.log(proportions)) / np.log(4))


def test_alpha_diversity_supports_aliases_and_multi_column_metrics() -> None:
    adata = ad.AnnData(
        X=np.array([[4, 3, 2, 1], [2, 0, 0, 0]], dtype=float),
        obs=pd.DataFrame(index=pd.Index(["s1", "s2"])),
        var=pd.DataFrame(index=pd.Index(["A", "B", "C", "D"])),
    )

    sobs = alpha_diversity(adata, "sobs")
    pielou = alpha_diversity(adata, "pielou")
    osd_result = alpha_diversity(adata, "osd")
    chao_ci = alpha_diversity(adata, "chao1-ci")

    assert sobs["sobs"].tolist() == [4.0, 1.0]
    assert pielou.columns.tolist() == ["sample_id", "pielou"]
    assert osd_result.columns.tolist() == ["sample_id", "observed_features", "singles", "doubles"]
    assert chao_ci.columns.tolist() == ["sample_id", "chao1_ci_lower", "chao1_ci_upper"]


def test_all_count_only_alpha_metrics_run_on_feature_table() -> None:
    adata = ad.AnnData(
        X=np.array([[4, 3, 2, 1], [2, 1, 1, 0]], dtype=float),
        obs=pd.DataFrame(index=pd.Index(["s1", "s2"])),
        var=pd.DataFrame(index=pd.Index(["A", "B", "C", "D"])),
    )
    count_only = set(available_alpha_metrics()) - {"faith_pd", "phydiv"}

    for metric in count_only:
        result = alpha_diversity(adata, metric)
        assert result.shape[0] == 2, metric
        assert result.columns[0] == "sample_id", metric


def test_phylogenetic_alpha_metrics_require_and_use_newick_tree() -> None:
    adata = ad.AnnData(
        X=np.array([[1, 0, 1], [1, 1, 0]], dtype=float),
        obs=pd.DataFrame(index=pd.Index(["s1", "s2"])),
        var=pd.DataFrame(index=pd.Index(["A", "B", "C"])),
    )
    tree = "((A:1,B:2):3,C:4)root;"

    faith = alpha_diversity(adata, "faith_pd", tree=tree)
    phydiv = alpha_diversity(adata, "phydiv", tree=tree)

    assert faith["faith_pd"].tolist() == [8.0, 6.0]
    assert phydiv["phydiv"].tolist() == [8.0, 6.0]

    with pytest.raises(MicrobiomeSuiteError, match="requires"):
        alpha_diversity(adata, "faith_pd")
