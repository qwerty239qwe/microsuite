from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from microsuite.api import (
    abundance_table,
    alpha_diversity,
    beta_diversity,
    normalize_table,
    pcoa,
    rarefy_table,
    read_table,
    shared_taxa_table,
    write_table,
)
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def test_python_sdk_facade_table_roundtrip_and_ecology(tmp_path: Path) -> None:
    table = tmp_path / "table.h5ad"
    write_table(fixture_adata(), table)
    adata = read_table(table)

    relative = normalize_table(adata, method="relative")
    abundance = abundance_table(adata, level="genus")
    shared = shared_taxa_table(adata, level="genus", group="body_site")
    rarefied = rarefy_table(adata, depth=10, seed=1)
    alpha = alpha_diversity(adata, metric="shannon")
    beta = beta_diversity(adata, metric="bray-curtis")
    coords = pcoa(beta, dimensions=2)

    assert np.allclose(np.asarray(relative.X).sum(axis=1), 1.0)
    assert isinstance(abundance, pd.DataFrame)
    assert "Lactobacillus" in abundance.columns
    assert shared["taxon"].tolist()
    assert np.asarray(rarefied.X).sum(axis=1).tolist() == [10.0, 10.0, 10.0, 10.0]
    assert "shannon" in alpha.columns
    assert beta.shape == (4, 4)
    assert coords.columns.tolist() == [
        "sample_id",
        "PC1",
        "PC2",
        "PC1_variance",
        "PC2_variance",
    ]
