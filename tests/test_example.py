from __future__ import annotations

from pathlib import Path

from microsuite.examples.moving_pictures import run_example


def test_moving_pictures_example(tmp_path: Path) -> None:
    run_example(tmp_path / "run")

    assert (tmp_path / "run" / "table.h5ad").exists()
    assert (tmp_path / "run" / "alpha-shannon.tsv").exists()
    assert (tmp_path / "run" / "beta-bray-curtis.tsv").exists()
    assert (tmp_path / "run" / "pcoa.tsv").exists()
    assert (tmp_path / "run" / "barplot-genus.png").exists()
    assert (tmp_path / "run" / "run.json").exists()
