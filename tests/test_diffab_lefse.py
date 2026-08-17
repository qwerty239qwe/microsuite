from __future__ import annotations

import json
import os
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.value_type import record_batch_correction
from microsuite.cli.app import app
from microsuite.diffab.lefse import run_lefse
from microsuite.io.h5ad import write_h5ad
from microsuite.methods.diff_abundance import diff_abundance


def _adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array(
            [
                [30, 3, 10],
                [28, 4, 11],
                [4, 30, 10],
                [3, 29, 12],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(
            {
                "group": ["A", "A", "B", "B"],
                "subclass": ["x", "y", "x", "y"],
            },
            index=["s1", "s2", "s3", "s4"],
        ),
        var=pd.DataFrame(index=["f1", "f2", "f3"]),
    )


def _fake_invoke(captured: dict[str, object], *, fail: bool = False):
    def fake(**kwargs: object) -> None:
        captured.update(kwargs)
        positional = kwargs["positional"]
        assert isinstance(positional, list)
        captured["params"] = json.loads(Path(positional[2]).read_text(encoding="utf-8"))
        if fail:
            raise MicrobiomeSuiteError("simulated failure")
        output = positional[-1]
        assert isinstance(output, Path)
        output.write_text("features\tscores\nf1\t3.5\n", encoding="utf-8")

    return fake


def test_lefse_forwards_design_reproducibility_and_thresholds(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}
    monkeypatch.setattr("microsuite.diffab.lefse.invoke_r_script", _fake_invoke(captured))
    output = tmp_path / "lefse.tsv"

    run_lefse(
        _adata(),
        output=output,
        group="group",
        subclass="subclass",
        reference="B",
        seed=77,
        kruskal_threshold=0.1,
        wilcoxon_threshold=0.2,
        lda_threshold=1.5,
        p_adjust_method="bh",
        trim_names=True,
    )

    assert captured["params"] == {
        "comparison": "A",
        "group": "group",
        "kruskal_threshold": 0.1,
        "lda_threshold": 1.5,
        "p_adjust_method": "BH",
        "reference": "B",
        "seed": 77,
        "subclass": "subclass",
        "trim_names": True,
        "wilcoxon_threshold": 0.2,
    }
    assert output.read_text(encoding="utf-8").startswith("features\tscores")
    assert json.loads((tmp_path / "lefse.tsv.params.json").read_text()) == captured["params"]


@pytest.mark.parametrize(
    ("mutator", "kwargs", "message"),
    [
        (lambda x: x, {"group": "missing"}, "Group column not found"),
        (lambda x: x, {"group": "group", "reference": "C"}, "reference"),
        (lambda x: x, {"group": "group", "subclass": "missing"}, "Subclass column"),
        (
            lambda x: x.obs.__setitem__("subclass", ["A1", "A2", "B1", "B2"]),
            {"group": "group", "subclass": "subclass"},
            "represented in both groups",
        ),
        (lambda x: x, {"group": "group", "seed": -1}, "seed"),
        (lambda x: x, {"group": "group", "kruskal_threshold": 1.1}, "kruskal"),
        (lambda x: x, {"group": "group", "wilcoxon_threshold": float("nan")}, "wilcoxon"),
        (lambda x: x, {"group": "group", "lda_threshold": -1}, "lda_threshold"),
        (lambda x: x, {"group": "group", "p_adjust_method": "bad"}, "p_adjust"),
        (
            lambda x: x.obs.__setitem__("group", ["", "A", "B", "B"]),
            {"group": "group"},
            "empty values",
        ),
        (
            lambda x: setattr(x, "X", np.array([[1, -1, 2]] * 4)),
            {"group": "group"},
            "non-negative",
        ),
        (
            lambda x: setattr(x, "X", np.array([[1, np.nan, 2]] * 4)),
            {"group": "group"},
            "non-finite",
        ),
        (
            lambda x: setattr(
                x, "X", np.array([[0, 0, 0], [1, 1, 1], [1, 1, 1], [1, 1, 1]])
            ),
            {"group": "group"},
            "zero total",
        ),
    ],
)
def test_lefse_rejects_invalid_inputs(
    tmp_path: Path, mutator, kwargs: dict[str, object], message: str
) -> None:
    adata = _adata()
    mutator(adata)
    with pytest.raises(MicrobiomeSuiteError, match=message):
        run_lefse(adata, output=tmp_path / "lefse.tsv", **kwargs)


def test_lefse_rejects_declared_clr(tmp_path: Path) -> None:
    adata = _adata()
    record_batch_correction(
        adata,
        value_type="clr",
        backend="plsda-batch",
        batch="batch",
        covariates=[],
        target=None,
    )
    with pytest.raises(MicrobiomeSuiteError, match="requires a table of type counts or relative"):
        run_lefse(adata, output=tmp_path / "lefse.tsv", group="group")


def test_lefse_force_is_transactional(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    output = tmp_path / "lefse.tsv"
    output.write_text("old\n", encoding="utf-8")
    monkeypatch.setattr(
        "microsuite.diffab.lefse.invoke_r_script", _fake_invoke({}, fail=True)
    )
    with pytest.raises(MicrobiomeSuiteError, match="simulated failure"):
        run_lefse(_adata(), output=output, group="group", force=True)
    assert output.read_text(encoding="utf-8") == "old\n"


def test_unified_method_dispatches_lefse(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    table = tmp_path / "table.h5ad"
    write_h5ad(_adata(), table)
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "microsuite.methods.diff_abundance.run_lefse",
        lambda adata, **kwargs: captured.update(kwargs),
    )
    diff_abundance(
        backend="lefse",
        table=table,
        group="group",
        output=tmp_path / "lefse.tsv",
        subclass="subclass",
        reference="B",
        seed=9,
    )
    assert captured["subclass"] == "subclass"
    assert captured["reference"] == "B"
    assert captured["seed"] == 9


def test_unified_method_rejects_cross_backend_options(tmp_path: Path) -> None:
    table = tmp_path / "table.h5ad"
    write_h5ad(_adata(), table)
    with pytest.raises(MicrobiomeSuiteError, match="cannot be used with --backend lefse"):
        diff_abundance(
            backend="lefse",
            table=table,
            group="group",
            formula="~ group",
            output=tmp_path / "lefse.tsv",
        )
    with pytest.raises(MicrobiomeSuiteError, match="cannot be used with --backend maaslin3"):
        diff_abundance(
            backend="maaslin3",
            table=table,
            formula="~ group",
            subclass="subclass",
            output=tmp_path / "maaslin3",
        )


def test_lefse_cli_accepts_hardened_options(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "microsuite.cli.method_stats_cmd.diff_abundance", lambda **kwargs: captured.update(kwargs)
    )
    result = CliRunner().invoke(
        app,
        [
            "diff_abundance", "--backend", "lefse", "--table", str(tmp_path / "x.h5ad"),
            "--group", "group", "--output", str(tmp_path / "lefse.tsv"),
            "--subclass", "subclass", "--reference", "B", "--seed", "55",
            "--kruskal-threshold", "0.1", "--wilcoxon-threshold", "0.2",
            "--lda-threshold", "1.5", "--p-adjust-method", "BY", "--trim-names",
        ],
    )
    assert result.exit_code == 0, result.output
    assert captured["subclass"] == "subclass"
    assert captured["reference"] == "B"
    assert captured["seed"] == 55
    assert captured["p_adjust_method"] == "BY"
    assert captured["trim_names"] is True


def test_lefse_r_script_parses() -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is unavailable")
    script = str(files("microsuite.diffab.r").joinpath("lefse.R"))
    result = subprocess.run(
        [rscript, "-e", f"invisible(parse(file={json.dumps(script)}))"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_LEFSE_SMOKE") != "1",
    reason="set MICROSUITE_RUN_LEFSE_SMOKE=1 for the live LEfSe test",
)
def test_lefse_live_recovers_planted_signal_and_is_reproducible(tmp_path: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is unavailable")
    package = subprocess.run(
        [rscript, "-e", 'quit(status=!requireNamespace("lefser", quietly=TRUE))'],
        capture_output=True,
    )
    if package.returncode != 0:
        pytest.skip("R package lefser is unavailable")

    rng = np.random.default_rng(812)
    n = 40
    group = np.array(["A"] * 20 + ["B"] * 20)
    matrix = rng.poisson(80, size=(n, 8)).astype(float)
    matrix[:, 0] = rng.poisson(np.where(group == "A", 15, 700))
    adata = ad.AnnData(
        X=matrix,
        obs=pd.DataFrame(
            {"group": group, "subclass": np.tile(np.repeat(["x", "y"], 10), 2)},
            index=[f"s{i}" for i in range(n)],
        ),
        var=pd.DataFrame(index=["planted_hit", *[f"background_{i}" for i in range(7)]]),
    )
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    run_lefse(adata, output=first, group="group", reference="A", seed=991)
    run_lefse(adata, output=second, group="group", reference="A", seed=991)
    blocked = tmp_path / "blocked.tsv"
    run_lefse(
        adata,
        output=blocked,
        group="group",
        subclass="subclass",
        reference="A",
        seed=991,
        p_adjust_method="BH",
    )

    result = pd.read_csv(first, sep="\t")
    assert "planted_hit" in set(result["features"])
    assert result.loc[result["features"] == "planted_hit", "scores"].iloc[0] > 0
    assert "planted_hit" in set(pd.read_csv(blocked, sep="\t")["features"])
    assert first.read_bytes() == second.read_bytes()
