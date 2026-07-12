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


def test_ordination_scatter_and_trajectory(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_braycurtis_ordination

    s = tmp_path / "scatter.png"
    plot_braycurtis_ordination(make_adata(), color_by="phase", output=s)  # scatter default
    assert s.exists() and s.stat().st_size > 0

    t = tmp_path / "traj.png"
    plot_braycurtis_ordination(  # numeric color + subject -> trajectory default
        make_adata(), color_by="time", subject="subject", output=t
    )
    assert t.exists() and t.stat().st_size > 0

    f = tmp_path / "facet.png"
    plot_braycurtis_ordination(
        make_adata(), color_by="phase", subject="subject", output=f, style="facet"
    )
    assert f.exists() and f.stat().st_size > 0


def test_ordination_trajectory_requires_subject(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_braycurtis_ordination

    with pytest.raises(MicrobiomeSuiteError, match="subject"):
        plot_braycurtis_ordination(
            make_adata(), color_by="phase", output=tmp_path / "x.png", style="trajectory"
        )


def test_cli_braycurtis_ordination(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(make_adata(), src)
    out = tmp_path / "ord.png"
    r = CliRunner().invoke(
        app,
        [
            "viz",
            "braycurtis-ordination",
            "--table",
            str(src),
            "--color-by",
            "time",
            "--subject",
            "subject",
            "-o",
            str(out),
        ],
    )
    assert r.exit_code == 0, r.stdout
    assert out.exists()


def test_plot_clr_by_group_styles(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_clr_by_group

    for style in ("boxplot", "heatmap", "violin"):
        out = tmp_path / f"clr_{style}.png"
        plot_clr_by_group(
            make_adata(), level="genus", group_by="phase", output=out, top_n=3, style=style
        )
        assert out.exists() and out.stat().st_size > 0


def test_plot_clr_by_group_bad_style(tmp_path: Path) -> None:
    from microsuite.viz.metadata import plot_clr_by_group

    with pytest.raises(MicrobiomeSuiteError, match="style"):
        plot_clr_by_group(
            make_adata(), level="genus", group_by="phase", output=tmp_path / "x.png", style="pie"
        )


def test_cli_clr_by_group_heatmap(tmp_path: Path) -> None:
    src = tmp_path / "t.h5ad"
    write_h5ad(make_adata(), src)
    out = tmp_path / "clr.png"
    r = CliRunner().invoke(
        app,
        [
            "viz",
            "clr-by-group",
            "--table",
            str(src),
            "--level",
            "genus",
            "--group-by",
            "phase",
            "--style",
            "heatmap",
            "-o",
            str(out),
        ],
    )
    assert r.exit_code == 0, r.stdout
    assert out.exists()
