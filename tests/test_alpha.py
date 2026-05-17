from __future__ import annotations

from pathlib import Path

import numpy as np

from microsuite.diversity.alpha import alpha_diversity
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
