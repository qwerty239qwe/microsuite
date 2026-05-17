from __future__ import annotations

from pathlib import Path

import numpy as np

from microsuite.diversity.beta import beta_diversity
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def test_beta_bray_curtis_known_value() -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    result = beta_diversity(adata, "bray-curtis")

    expected = 18 / 22
    assert np.isclose(result.loc["L1S8", "L1S57"], expected)
    assert np.isclose(result.loc["L1S57", "L1S8"], expected)
    assert np.isclose(result.loc["L1S8", "L1S8"], 0.0)


def test_beta_jaccard_known_value() -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    result = beta_diversity(adata, "jaccard")

    assert np.isclose(result.loc["L1S8", "L1S57"], 0.75)
