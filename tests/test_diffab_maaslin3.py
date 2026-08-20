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
from microsuite.cli.app import app
from microsuite.diffab.maaslin3 import run_maaslin3
from microsuite.io.h5ad import write_h5ad
from microsuite.methods.diff_abundance import diff_abundance


def _adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array([[10, 0, 20], [12, 1, 19], [30, 8, 18], [32, 9, 21]]),
        obs=pd.DataFrame(
            {
                "group": ["A", "A", "B", "B"],
                "batch": ["x", "x", "y", "y"],
                "subject": ["p1", "p2", "p1", "p2"],
            },
            index=["s1", "s2", "s3", "s4"],
        ),
        var=pd.DataFrame(index=["f1", "f2", "f3"]),
    )


def _fake_invoke(captured: dict[str, object]):
    def fake(**kwargs: object) -> None:
        captured.update(kwargs)
        positional = kwargs["positional"]
        assert isinstance(positional, list)
        captured["params"] = json.loads(Path(positional[2]).read_text(encoding="utf-8"))
        output = positional[-1]
        assert isinstance(output, Path)
        output.mkdir()
        (output / "abundance_results.tsv").write_text(
            "feature\tmodel\nf1\tabundance\n", encoding="utf-8"
        )
        (output / "prevalence_results.tsv").write_text(
            "feature\tmodel\nf2\tprevalence\n", encoding="utf-8"
        )
        (output / "all_results.tsv").write_text(
            "feature\tmodel\nf1\tabundance\nf2\tprevalence\n", encoding="utf-8"
        )

    return fake


def test_maaslin3_forwards_full_formula_and_options(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}
    monkeypatch.setattr("microsuite.diffab.maaslin3.invoke_r_script", _fake_invoke(captured))

    output = tmp_path / "results"
    run_maaslin3(
        _adata(),
        output=output,
        formula="~ batch + group + (1 | subject)",
        normalization="clr",
        transform="none",
        min_prevalence=0.25,
        min_abundance=0.01,
    )

    params = captured["params"]
    assert params == {
        "formula": "~ batch + group + (1 | subject)",
        "normalization": "CLR",
        "transform": "NONE",
        "min_prevalence": 0.25,
        "min_abundance": 0.01,
    }
    assert (output / "abundance_results.tsv").is_file()
    assert (output / "prevalence_results.tsv").is_file()


def test_maaslin3_forwards_reference_levels(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A named reference reaches params; omitting it leaves params unchanged."""
    captured: dict[str, object] = {}
    monkeypatch.setattr("microsuite.diffab.maaslin3.invoke_r_script", _fake_invoke(captured))

    run_maaslin3(
        _adata(),
        output=tmp_path / "with-reference",
        formula="~ batch + group",
        reference=" group,B ; batch,x ",
    )
    assert captured["params"]["reference"] == "group,B;batch,x"

    run_maaslin3(
        _adata(),
        output=tmp_path / "without-reference",
        formula="~ batch + group",
    )
    assert "reference" not in captured["params"]


@pytest.mark.parametrize(
    ("reference", "message"),
    [
        ("", "cannot be empty"),
        ("group", "expected 'column,level'"),
        ("group,A;group,B", "specified more than once"),
        ("missing,A", "column not found"),
    ],
)
def test_maaslin3_rejects_invalid_reference_contract(
    tmp_path: Path, reference: str, message: str
) -> None:
    with pytest.raises(MicrobiomeSuiteError, match=message):
        run_maaslin3(
            _adata(),
            output=tmp_path / "results",
            formula="~ group",
            reference=reference,
        )


def test_maaslin3_combines_fixed_and_random_formulas(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}
    monkeypatch.setattr("microsuite.diffab.maaslin3.invoke_r_script", _fake_invoke(captured))

    run_maaslin3(
        _adata(),
        output=tmp_path / "results",
        fix_formula="batch + group",
        rand_formula="(1 | subject)",
    )

    params = captured["params"]
    assert params["formula"] == "~ batch + group + (1 | subject)"


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({}, "Provide --formula"),
        ({"group": "group", "fix_formula": "batch"}, "either --group or --fix-formula"),
        ({"formula": "~ group", "rand_formula": "(1|subject)"}, "cannot be combined"),
        ({"group": "missing"}, "Group column not found"),
        ({"group": "group", "normalization": "bad"}, "normalization"),
        ({"group": "group", "transform": "bad"}, "transform"),
        ({"group": "group", "min_prevalence": 1.1}, "min_prevalence"),
        ({"group": "group", "min_prevalence": float("inf")}, "min_prevalence"),
        ({"group": "group", "min_abundance": -1}, "min_abundance"),
        ({"group": "group", "min_abundance": float("nan")}, "min_abundance"),
    ],
)
def test_maaslin3_rejects_invalid_designs_and_options(
    tmp_path: Path, kwargs: dict[str, object], message: str
) -> None:
    with pytest.raises(MicrobiomeSuiteError, match=message):
        run_maaslin3(_adata(), output=tmp_path / "results", **kwargs)


def test_maaslin3_force_replaces_output_only_after_success(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    output = tmp_path / "results"
    output.mkdir()
    (output / "old.txt").write_text("old", encoding="utf-8")
    monkeypatch.setattr("microsuite.diffab.maaslin3.invoke_r_script", _fake_invoke({}))

    run_maaslin3(_adata(), output=output, group="group", force=True)

    assert not (output / "old.txt").exists()
    assert (output / "abundance_results.tsv").is_file()


def test_unified_method_dispatches_maaslin3(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.h5ad"
    write_h5ad(_adata(), table)
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "microsuite.methods.diff_abundance.run_maaslin3",
        lambda adata, **kwargs: captured.update(kwargs),
    )

    diff_abundance(
        backend="maaslin3",
        table=table,
        group=None,
        formula="~ batch + group + (1 | subject)",
        output=tmp_path / "results",
        normalization="CLR",
        reference="group,B;batch,x",
    )

    assert captured["formula"] == "~ batch + group + (1 | subject)"
    assert captured["normalization"] == "CLR"
    assert captured["reference"] == "group,B;batch,x"
    assert captured["force"] is False


def test_maaslin3_cli_accepts_formula_without_group(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.h5ad"
    write_h5ad(_adata(), table)
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "microsuite.cli.method_stats_cmd.diff_abundance",
        lambda **kwargs: captured.update(kwargs),
    )

    result = CliRunner().invoke(
        app,
        [
            "diff_abundance",
            "--backend",
            "maaslin3",
            "--table",
            str(table),
            "--output",
            str(tmp_path / "results"),
            "--formula",
            "~ batch + group + (1 | subject)",
            "--min-prevalence",
            "0.2",
            "--reference",
            "group,B;batch,x",
        ],
    )

    assert result.exit_code == 0, result.output
    assert captured["group"] is None
    assert captured["formula"] == "~ batch + group + (1 | subject)"
    assert captured["min_prevalence"] == 0.2
    assert captured["reference"] == "group,B;batch,x"


def test_maaslin3_r_script_parses() -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is unavailable")
    script = str(files("microsuite.diffab.r").joinpath("maaslin3.R"))
    result = subprocess.run(
        [rscript, "-e", f"invisible(parse(file={json.dumps(script)}))"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_MAASLIN3_SMOKE") != "1",
    reason="set MICROSUITE_RUN_MAASLIN3_SMOKE=1 for the live MaAsLin 3 test",
)
def test_maaslin3_live_recovers_planted_abundance_and_prevalence(
    tmp_path: Path,
) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is unavailable")
    package = subprocess.run(
        [rscript, "-e", 'quit(status=!requireNamespace("maaslin3", quietly=TRUE))'],
        capture_output=True,
    )
    if package.returncode != 0:
        pytest.skip("R package maaslin3 is unavailable")

    n = 80
    group = np.array(["A"] * 40 + ["B"] * 40)
    subject = np.repeat([f"p{i}" for i in range(20)], 4)
    rng = np.random.default_rng(404)
    values = rng.poisson(100, size=(n, 10)).astype(float)
    values[:, 0] = rng.poisson(np.where(group == "A", 15, 500))
    values[:, 1] = 0
    values[:4, 1] = rng.poisson(100, size=4)
    values[40:76, 1] = rng.poisson(100, size=36)
    adata = ad.AnnData(
        X=values,
        obs=pd.DataFrame({"group": group, "subject": subject}, index=[f"s{i}" for i in range(n)]),
        var=pd.DataFrame(index=["abundance_hit", "prevalence_hit", *[f"bg{i}" for i in range(8)]]),
    )

    output = tmp_path / "results"
    run_maaslin3(adata, output=output, formula="~ group + (1 | subject)")
    abundance = pd.read_csv(output / "abundance_results.tsv", sep="\t")
    prevalence = pd.read_csv(output / "prevalence_results.tsv", sep="\t")

    assert set(abundance["model"]) == {"abundance"}
    assert set(prevalence["model"]) == {"prevalence"}
    abundance_hit = abundance[abundance["feature"] == "abundance_hit"]
    prevalence_hit = prevalence[prevalence["feature"] == "prevalence_hit"]
    assert abundance_hit["qval_individual"].min() < 0.05
    assert prevalence_hit["qval_individual"].min() < 0.05


@pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_MAASLIN3_SMOKE") != "1",
    reason="set MICROSUITE_RUN_MAASLIN3_SMOKE=1 for the live MaAsLin 3 test",
)
def test_maaslin3_live_multilevel_reference_controls_baseline(tmp_path: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is unavailable")
    package = subprocess.run(
        [rscript, "-e", 'quit(status=!requireNamespace("maaslin3", quietly=TRUE))'],
        capture_output=True,
    )
    if package.returncode != 0:
        pytest.skip("R package maaslin3 is unavailable")

    group = np.repeat(["A", "B", "C"], 20)
    rng = np.random.default_rng(405)
    values = rng.poisson(100, size=(60, 4)).astype(float)
    values[:, 0] += np.select([group == "A", group == "B"], [120, 60], default=0)
    adata = ad.AnnData(
        X=values,
        obs=pd.DataFrame(
            {"group": group, "unused_subject": [f"p{i}" for i in range(60)]},
            index=[f"s{i}" for i in range(60)],
        ),
        var=pd.DataFrame(index=["hit", "bg1", "bg2", "bg3"]),
    )

    explicit = tmp_path / "explicit"
    explicit_run = tmp_path / "explicit-run"
    run_maaslin3(
        adata,
        output=explicit,
        formula="~ group",
        reference="group,C",
        run_dir=explicit_run,
    )
    explicit_results = pd.read_csv(explicit / "all_results.tsv", sep="\t")
    group_values = set(explicit_results.loc[explicit_results["metadata"] == "group", "value"])
    assert group_values == {"A", "B"}
    stderr = (explicit_run / "stderr.log").read_text(encoding="utf-8")
    assert "'group' reference level = 'C'" in stderr
    assert "unused_subject" not in stderr

    default = tmp_path / "default"
    run_maaslin3(adata, output=default, formula="~ group")
    default_results = pd.read_csv(default / "all_results.tsv", sep="\t")
    default_values = set(default_results.loc[default_results["metadata"] == "group", "value"])
    assert default_values == {"B", "C"}

    counts_path = tmp_path / "counts.tsv"
    metadata_path = tmp_path / "metadata.tsv"
    params_path = tmp_path / "bad-params.json"
    pd.DataFrame(values.T, index=adata.var_names, columns=adata.obs_names).to_csv(
        counts_path, sep="\t"
    )
    pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")
    params_path.write_text(
        json.dumps(
            {
                "formula": "~ group",
                "reference": "group",
                "normalization": "TSS",
                "transform": "LOG",
                "min_prevalence": 0,
                "min_abundance": 0,
            }
        ),
        encoding="utf-8",
    )
    script = str(files("microsuite.diffab.r").joinpath("maaslin3.R"))
    malformed = subprocess.run(
        [rscript, script, counts_path, metadata_path, params_path, tmp_path / "bad-output"],
        capture_output=True,
        text=True,
    )
    assert malformed.returncode != 0
    assert "expected 'column,level' pairs" in malformed.stderr
