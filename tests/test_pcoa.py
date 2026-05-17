from __future__ import annotations

from pathlib import Path

from microsuite.diversity.beta import beta_diversity
from microsuite.io.tsv import read_tsv
from microsuite.ordination.pcoa import pcoa

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def test_pcoa_outputs_sample_ids_and_dimensions() -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    result = pcoa(beta_diversity(adata, "bray-curtis"), dimensions=2)

    assert result.columns.tolist() == [
        "sample_id",
        "PC1",
        "PC2",
        "PC1_variance",
        "PC2_variance",
    ]
    assert result["sample_id"].tolist() == ["L1S8", "L1S57", "L1S76", "L2S155"]
