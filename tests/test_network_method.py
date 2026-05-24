from __future__ import annotations

import json
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_tsv
from microsuite.methods.network import correlation_network, network, sparcc_network

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def fixture_table(tmp_path: Path) -> Path:
    table = tmp_path / "table.h5ad"
    write_h5ad(fixture_adata(), table)
    return table


def test_network_runner_scripts_are_external_assets() -> None:
    r_script = ROOT / "scripts" / "r" / "spieceasi_network.R"
    packaged_r = files("microsuite.networks.r").joinpath("spieceasi_network.R")
    packaged_julia = files("microsuite.networks.julia").joinpath("flashweave_network.jl")

    assert r_script.exists()
    assert "SpiecEasi::spiec.easi" in r_script.read_text(encoding="utf-8")
    assert packaged_r.is_file()
    assert packaged_r.read_text(encoding="utf-8") == r_script.read_text(encoding="utf-8")
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


def test_sparcc_network_uses_clr_backend_label() -> None:
    edges = sparcc_network(fixture_adata(), min_abs_weight=0.0, min_prevalence=0.0)

    assert not edges.empty
    assert edges["backend"].unique().tolist() == ["sparcc"]
    assert edges["p_value"].isna().all()


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
