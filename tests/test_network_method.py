from __future__ import annotations

import json
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

import microsuite.methods.network as network_method
from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_tsv
from microsuite.methods._sparcc import SparCCResult
from microsuite.methods.network import correlation_network, network, sparcc_network

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def fixture_table(tmp_path: Path) -> Path:
    table = tmp_path / "table.h5ad"
    write_h5ad(fixture_adata(), table)
    return table


def test_network_runner_scripts_are_external_assets() -> None:
    packaged_r = files("microsuite.networks.r").joinpath("spieceasi_network.R")
    packaged_julia = files("microsuite.networks.julia").joinpath("flashweave_network.jl")

    assert packaged_r.is_file()
    assert "SpiecEasi::spiec.easi" in packaged_r.read_text(encoding="utf-8")
    assert packaged_julia.is_file()
    assert "learn_network" in packaged_julia.read_text(encoding="utf-8")


def test_native_correlation_network_outputs_edges() -> None:
    edges = correlation_network(
        fixture_adata(),
        method="pearson",
        transform="counts",
        min_abs_weight=0.0,
        min_prevalence=0.0,
    )

    assert {"source", "target", "weight", "p_value", "backend"}.issubset(edges.columns)
    assert edges["backend"].unique().tolist() == ["native-correlation"]
    assert edges["abs_weight"].is_monotonic_decreasing


def test_sparcc_network_filters_raw_counts_and_uses_estimated_correlations(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    adata = ad.AnnData(
        np.array(
            [
                [5, 1, 0, 2],
                [0, 0, 4, 1],
                [3, 0, 2, 0],
                [1, 0, 1, 5],
            ]
        ),
        var=pd.DataFrame(index=["zeta", "filtered", "alpha", "mu"]),
    )
    correlation = np.array(
        [
            [1.0, 0.5, -0.8],
            [0.5, 1.0, 0.7],
            [-0.8, 0.7, 1.0],
        ]
    )
    calls: list[tuple[np.ndarray, dict[str, object]]] = []

    def fake_estimate(counts: np.ndarray, **parameters: object) -> SparCCResult:
        calls.append((counts.copy(), parameters))
        return SparCCResult(covariance=np.eye(3), correlation=correlation)

    def forbidden_call(*args: object, **kwargs: object) -> None:
        raise AssertionError("SparCC must not use CLR transformation or Pearson correlation")

    monkeypatch.setattr(network_method, "estimate_sparcc", fake_estimate)
    monkeypatch.setattr(network_method, "_transform_counts", forbidden_call)
    monkeypatch.setattr(network_method.stats, "pearsonr", forbidden_call)

    edges = sparcc_network(
        adata,
        min_abs_weight=0.5,
        min_prevalence=0.5,
        pseudocount=0.25,
        iterations=7,
        inner_iterations=4,
        exclusion_threshold=0.2,
        seed=99,
    )
    limited = sparcc_network(
        adata,
        min_abs_weight=0.5,
        min_prevalence=0.5,
        top_n=2,
        pseudocount=0.25,
        iterations=7,
        inner_iterations=4,
        exclusion_threshold=0.2,
        seed=99,
    )

    expected_counts = np.array([[5, 0, 2], [0, 4, 1], [3, 2, 0], [1, 1, 5]])
    for counts, parameters in calls:
        np.testing.assert_array_equal(counts, expected_counts)
        assert parameters == {
            "iterations": 7,
            "inner_iterations": 4,
            "exclusion_threshold": 0.2,
            "pseudocount": 0.25,
            "seed": 99,
        }
    assert edges[["source", "target"]].values.tolist() == [
        ["zeta", "mu"],
        ["alpha", "mu"],
        ["zeta", "alpha"],
    ]
    np.testing.assert_array_equal(edges["weight"].to_numpy(), [-0.8, 0.7, 0.5])
    np.testing.assert_array_equal(edges["abs_weight"].to_numpy(), [0.8, 0.7, 0.5])
    assert edges.columns.tolist() == [
        "source",
        "target",
        "weight",
        "abs_weight",
        "p_value",
        "method",
        "backend",
    ]
    assert edges["p_value"].isna().all()
    assert edges["method"].unique().tolist() == ["sparcc"]
    assert edges["backend"].unique().tolist() == ["sparcc"]
    pd.testing.assert_frame_equal(limited, edges.head(2).reset_index(drop=True))


def test_sparcc_network_empty_result_has_exact_edge_schema(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    adata = ad.AnnData(np.array([[2, 1, 3], [1, 4, 2]]))
    correlation = np.array(
        [
            [1.0, 0.1, -0.1],
            [0.1, 1.0, 0.1],
            [-0.1, 0.1, 1.0],
        ]
    )
    monkeypatch.setattr(
        network_method,
        "estimate_sparcc",
        lambda *args, **kwargs: SparCCResult(
            covariance=np.eye(3),
            correlation=correlation,
        ),
    )

    edges = sparcc_network(adata, min_abs_weight=0.2, min_prevalence=0.0)

    assert edges.empty
    assert edges.columns.tolist() == [
        "source",
        "target",
        "weight",
        "abs_weight",
        "p_value",
        "method",
        "backend",
    ]
    assert edges.dtypes.astype(str).tolist() == [
        "string",
        "string",
        "float64",
        "float64",
        "float64",
        "string",
        "string",
    ]


@pytest.mark.parametrize(
    "matrix",
    [
        np.array([[0.5, 0.25, 0.25], [0.2, 0.3, 0.5]]),
        np.array([[1.5, 2.0, 3.0], [2.0, 1.25, 4.0]]),
    ],
)
def test_sparcc_network_rejects_normalized_or_fractional_anndata(matrix: np.ndarray) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="raw counts"):
        sparcc_network(ad.AnnData(matrix), min_abs_weight=0.0, min_prevalence=0.0)


def test_sparcc_network_is_deterministic_for_the_same_seed() -> None:
    adata = ad.AnnData(np.array([[4, 0, 2], [1, 3, 1], [0, 2, 5], [2, 1, 2]]))

    first = sparcc_network(
        adata,
        min_abs_weight=0.0,
        min_prevalence=0.0,
        iterations=3,
        inner_iterations=2,
        seed=12,
    )
    repeated = sparcc_network(
        adata,
        min_abs_weight=0.0,
        min_prevalence=0.0,
        iterations=3,
        inner_iterations=2,
        seed=12,
    )

    pd.testing.assert_frame_equal(first, repeated)


def test_network_sparcc_forwards_tuning_parameters(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_sparcc(adata: ad.AnnData, **parameters: object) -> pd.DataFrame:
        captured.update(parameters)
        return pd.DataFrame({"source": pd.Series(dtype="string")})

    monkeypatch.setattr(network_method, "sparcc_network", fake_sparcc)

    network(
        backend="sparcc",
        table=fixture_table(tmp_path),
        output=tmp_path / "sparcc.tsv",
        pseudocount=0.75,
        iterations=9,
        inner_iterations=6,
        exclusion_threshold=0.35,
        seed=41,
    )

    assert captured == {
        "min_abs_weight": 0.3,
        "min_prevalence": 0.1,
        "top_n": None,
        "pseudocount": 0.75,
        "iterations": 9,
        "inner_iterations": 6,
        "exclusion_threshold": 0.35,
        "seed": 41,
    }


def test_network_cli_native_correlation_writes_edge_list(tmp_path: Path) -> None:
    table = fixture_table(tmp_path)
    output = tmp_path / "network.tsv"

    result = CliRunner().invoke(
        app,
        [
            "network",
            "infer",
            "--backend",
            "native-correlation",
            "--table",
            str(table),
            "-o",
            str(output),
            "--min-abs-weight",
            "0",
        ],
    )

    assert result.exit_code == 0, result.stdout
    edges = pd.read_csv(output, sep="\t")
    assert "weight" in edges.columns


def test_spieceasi_builds_command_and_logs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    table = fixture_table(tmp_path)
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript" if name == "Rscript" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    run_dir = tmp_path / "run"
    network(
        backend="spieceasi",
        table=table,
        output=tmp_path / "spieceasi.tsv",
        spieceasi_method="glasso",
        lambda_min_ratio=0.05,
        nlambda=5,
        run_dir=run_dir,
    )

    assert calls
    command = calls[0]
    assert command[0] == "Rscript"
    assert command[1].endswith("microsuite/networks/r/spieceasi_network.R") or command[1].endswith(
        "microsuite\\networks\\r\\spieceasi_network.R"
    )
    assert command[-3:] == ["glasso", "0.05", "5"]
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "network"
    assert run["backend"] == "spieceasi"


def test_flashweave_builds_command(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    table = fixture_table(tmp_path)
    calls: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "julia" if name == "julia" else None)
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: (
            calls.append(command) or subprocess.CompletedProcess(command, 0, "ok\n", "")
        ),
    )

    network(
        backend="flashweave",
        table=table,
        output=tmp_path / "flashweave.edgelist",
        sensitive=True,
        heterogeneous=True,
    )

    assert calls
    command = calls[0]
    assert command[0] == "julia"
    assert command[1].endswith("microsuite/networks/julia/flashweave_network.jl") or command[
        1
    ].endswith("microsuite\\networks\\julia\\flashweave_network.jl")
    assert command[-3:] == [str(tmp_path / "flashweave.edgelist"), "true", "true"]


@pytest.mark.parametrize(
    ("backend", "message"),
    [("spieceasi", "Rscript"), ("flashweave", "Julia")],
)
def test_network_external_backends_report_missing_runtime(
    backend: str, message: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(shutil, "which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match=message):
        network(
            backend=backend,
            table=fixture_table(tmp_path),
            output=tmp_path / "network.tsv",
        )
