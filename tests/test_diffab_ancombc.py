from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diffab.ancombc import run_ancombc


def _adata() -> ad.AnnData:
    rng = np.random.default_rng(0)
    X = rng.integers(0, 40, size=(6, 4)).astype(float)
    obs = pd.DataFrame(
        {
            "phase": ["pre", "pre", "post", "post", "pre", "post"],
            "subject": ["s1", "s2", "s1", "s2", "s3", "s3"],
        },
        index=[f"S{i}" for i in range(6)],
    )
    return ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=[f"F{i}" for i in range(4)]))


def _capture(monkeypatch) -> dict:
    captured: dict = {}

    def fake_run_command(command, **kw):
        captured["command"] = command
        captured["params"] = json.loads(Path(command[4]).read_text())
        # counts/metadata written to command[2]/command[3]
        captured["counts_exists"] = Path(command[2]).exists()

    monkeypatch.setattr("microsuite.runtime.r_backend.run_command", fake_run_command)
    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/Rscript")
    return captured


def test_run_ancombc_params_defaults(tmp_path, monkeypatch) -> None:
    captured = _capture(monkeypatch)
    run_ancombc(_adata(), output=tmp_path / "out.tsv", group="phase")
    p = captured["params"]
    assert p["fix_formula"] == "phase"  # defaulted from --group
    assert p["group"] == "phase"
    assert p["rand_formula"] is None
    # ANCOM-BC2 native defaults, not force-all-off
    assert p["prv_cut"] == 0.10 and p["lib_cut"] == 0
    assert p["pseudo_sens"] is True and p["p_adj_method"] == "BH"
    assert p["global"] is False and p["n_cl"] == 1
    assert captured["command"][1].endswith("ancombc.R")
    assert captured["counts_exists"]


def test_run_ancombc_fix_formula_overrides_group_and_rand(tmp_path, monkeypatch) -> None:
    captured = _capture(monkeypatch)
    run_ancombc(
        _adata(),
        output=tmp_path / "out.tsv",
        group="phase",
        fix_formula="phase*subject",
        rand_formula="(1|subject)",
        prv_cut=0.2,
        struc_zero=True,
        global_test=True,
        n_cl=4,
        reference={"phase": "pre"},
    )
    p = captured["params"]
    assert p["fix_formula"] == "phase*subject"  # overrides group
    assert p["rand_formula"] == "(1|subject)"
    assert p["reference"] == {"phase": "pre"}
    assert p["prv_cut"] == 0.2 and p["struc_zero"] is True
    assert p["global"] is True and p["n_cl"] == 4


def test_run_ancombc_requires_a_formula(tmp_path, monkeypatch) -> None:
    _capture(monkeypatch)
    with pytest.raises(MicrobiomeSuiteError, match="fix-formula|group"):
        run_ancombc(_adata(), output=tmp_path / "out.tsv")


def test_run_ancombc_missing_reference_column_raises(tmp_path, monkeypatch) -> None:
    _capture(monkeypatch)
    with pytest.raises(MicrobiomeSuiteError, match="nope"):
        run_ancombc(_adata(), output=tmp_path / "out.tsv", group="phase", reference={"nope": "x"})


def test_cli_ancombc_threads_options(tmp_path, monkeypatch) -> None:
    from typer.testing import CliRunner

    from microsuite.cli.app import app
    from microsuite.io.h5ad import write_h5ad

    src = tmp_path / "t.h5ad"
    write_h5ad(_adata(), src)
    captured: dict = {}

    def fake_run_ancombc(adata, **kwargs):
        captured.update(kwargs)

    monkeypatch.setattr("microsuite.cli.diffab_cmd.run_ancombc", fake_run_ancombc)
    r = CliRunner().invoke(
        app,
        [
            "diffab",
            "ancombc",
            str(src),
            "--group",
            "phase",
            "--rand-formula",
            "(1|subject)",
            "--reference",
            "phase=pre",
            "--prv-cut",
            "0.2",
            "--global",
            "--n-cl",
            "4",
            "-o",
            str(tmp_path / "o.tsv"),
        ],
    )
    assert r.exit_code == 0, r.stdout
    assert captured["group"] == "phase"
    assert captured["rand_formula"] == "(1|subject)"
    assert captured["reference"] == {"phase": "pre"}
    assert captured["prv_cut"] == 0.2
    assert captured["global_test"] is True and captured["n_cl"] == 4


def test_cli_ancombc_runtime_docker_image_reaches_backend(tmp_path, monkeypatch) -> None:
    from typer.testing import CliRunner

    from microsuite.cli.app import app
    from microsuite.io.h5ad import write_h5ad

    src = tmp_path / "t.h5ad"
    write_h5ad(_adata(), src)
    captured: dict = {}
    monkeypatch.setattr(
        "microsuite.cli.diffab_cmd.run_ancombc",
        lambda adata, **kwargs: captured.update(kwargs),
    )
    r = CliRunner().invoke(
        app,
        [
            "diffab",
            "ancombc",
            str(src),
            "--group",
            "phase",
            "--runtime",
            "docker",
            "--image",
            "my/img:1",
            "-o",
            str(tmp_path / "o.tsv"),
        ],
    )
    assert r.exit_code == 0, r.stdout
    assert captured["runtime"] == "docker" and captured["image"] == "my/img:1"
