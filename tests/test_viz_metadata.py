from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.viz.metadata import plot_taxa_by_group


def make_adata() -> ad.AnnData:
    # 6 samples, 4 features across 2 genera; obs has phase/time/subject
    rng = np.random.default_rng(0)
    X = rng.integers(1, 50, size=(6, 4)).astype(float)
    var = pd.DataFrame(
        {"genus": ["Bacteroides", "Bacteroides", "Prevotella", "Faecalibacterium"]},
        index=["F1", "F2", "F3", "F4"],
    )
    obs = pd.DataFrame(
        {
            "phase": ["pre", "pre", "post", "post", "pre", "post"],
            "time": [0, 0, 7, 7, 0, 7],
            "subject": ["s1", "s2", "s1", "s2", "s3", "s3"],
        },
        index=[f"S{i}" for i in range(6)],
    )
    return ad.AnnData(X=X, obs=obs, var=var)


def test_plot_taxa_by_group_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "taxa.png"
    plot_taxa_by_group(make_adata(), level="genus", group_by="phase", output=out, top_n=2)
    assert out.exists() and out.stat().st_size > 0


def test_plot_taxa_by_group_bad_column(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="nope"):
        plot_taxa_by_group(make_adata(), level="genus", group_by="nope", output=tmp_path / "x.png")


def test_plot_taxa_by_group_bad_level(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        plot_taxa_by_group(
            make_adata(), level="species", group_by="phase", output=tmp_path / "x.png"
        )


def test_cli_taxa_by_group(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(make_adata(), src)
    out = tmp_path / "taxa.png"
    r = CliRunner().invoke(
        app,
        [
            "viz",
            "taxa-by-group",
            "--table",
            str(src),
            "--level",
            "genus",
            "--group-by",
            "phase",
            "-o",
            str(out),
        ],
    )
    assert r.exit_code == 0, r.stdout
    assert out.exists()
