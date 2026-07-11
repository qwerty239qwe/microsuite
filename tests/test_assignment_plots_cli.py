from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from typer.testing import CliRunner

from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.viz.assignment import (
    plot_assigned_asv_by_rank,
    plot_assigned_reads_by_rank,
    plot_deepest_rank,
)


def _fixture() -> ad.AnnData:
    var = pd.DataFrame(
        {
            "kingdom": ["Bacteria", "Bacteria", ""],
            "phylum": ["Firmicutes", "Bacteroidetes", ""],
            "class": ["Bacilli", "Bacteroidia", ""],
            "order": ["Lactobacillales", "Bacteroidales", ""],
            "family": ["Lactobacillaceae", "Prevotellaceae", ""],
            "genus": ["Lactobacillus", "Prevotella", ""],
            "species": ["L. casei", "", ""],
        },
        index=["F1", "F2", "F3"],
    )
    X = np.array([[10.0, 5.0, 1.0], [0.0, 3.0, 2.0]])
    return ad.AnnData(X=X, obs=pd.DataFrame(index=["s1", "s2"]), var=var)


def test_plot_functions_write_pngs(tmp_path: Path) -> None:
    adata = _fixture()
    for fn, name in (
        (plot_assigned_asv_by_rank, "asv.png"),
        (plot_assigned_reads_by_rank, "reads.png"),
        (plot_deepest_rank, "deepest.png"),
    ):
        out = tmp_path / name
        fn(adata, out)
        assert out.exists() and out.stat().st_size > 0


def test_plot_reads_single_sample(tmp_path: Path) -> None:
    adata = _fixture()[[0]].copy()  # one sample
    out = tmp_path / "reads.png"
    plot_assigned_reads_by_rank(adata, out)
    assert out.exists() and out.stat().st_size > 0


def test_cli_assignment_summary_and_plots(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(_fixture(), src)
    runner = CliRunner()

    summary = tmp_path / "summary.tsv"
    r1 = runner.invoke(app, ["tax_assignment_summary", "--table", str(src), "-o", str(summary)])
    assert r1.exit_code == 0, r1.stdout
    assert summary.exists()

    plots = tmp_path / "plots"
    r2 = runner.invoke(
        app, ["tax_assignment_plots", "--table", str(src), "--output-dir", str(plots)]
    )
    assert r2.exit_code == 0, r2.stdout
    for name in ("assigned_asv_by_rank.png", "assigned_reads_by_rank.png", "deepest_rank.png"):
        assert (plots / name).exists()
