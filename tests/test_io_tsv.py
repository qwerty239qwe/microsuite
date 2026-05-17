from __future__ import annotations

from pathlib import Path

import pandas as pd

from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def test_tsv_import_ann_data_shape_and_taxonomy() -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")

    assert adata.shape == (4, 4)
    assert list(adata.obs_names) == ["L1S8", "L1S57", "L1S76", "L2S155"]
    assert list(adata.var_names) == ["f1", "f2", "f3", "f4"]
    assert pd.DataFrame(adata.var).loc["f1", "genus"] == "Lactobacillus"
    assert adata.uns["microsuite"]["importer"] == "tsv"
