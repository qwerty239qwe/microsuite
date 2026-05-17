from __future__ import annotations

from pathlib import Path

from microsuite.io.tsv import read_tsv
from microsuite.viz.barplot import taxonomy_barplot

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def test_barplot_writes_png(tmp_path: Path) -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    output = tmp_path / "barplot.png"

    taxonomy_barplot(adata, level="genus", output=output)

    assert output.exists()
    assert output.stat().st_size > 0
