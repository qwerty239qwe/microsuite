from __future__ import annotations

import json
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path
from types import SimpleNamespace

import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_tsv
from microsuite.methods.diff_abundance import diff_abundance

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_table(tmp_path: Path) -> Path:
    table = tmp_path / "table.h5ad"
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    write_h5ad(adata, table)
    return table


@pytest.mark.parametrize("script_name", ["ancombc", "aldex2", "maaslin2", "maaslin3", "lefse"])
def test_r_diffab_scripts_are_external_assets(script_name: str) -> None:
    packaged_script = files("microsuite.diffab.r").joinpath(f"{script_name}.R")

    assert packaged_script.is_file()
    assert "commandArgs" in packaged_script.read_text(encoding="utf-8")


def test_diff_abundance_ancombc_missing_rscript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="Rscript"):
        diff_abundance(
            backend="ancombc",
            table=table,
            group="treatment",
            output=tmp_path / "diff.tsv",
        )


def test_ancombc_invokes_external_script_path(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    commands: list[list[str]] = []
    captured_params: list[dict[str, object]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript")

    def fake_run(command: list[str], **kwargs: object) -> object:
        commands.append(command)
        captured_params.append(json.loads(Path(command[4]).read_text(encoding="utf-8")))
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)

    diff_abundance(
        backend="ancombc",
        table=table,
        group="treatment",
        output=tmp_path / "diff.tsv",
    )

    assert commands
    command = commands[0]
    assert command[0] == "Rscript"
    assert command[1].endswith("microsuite/diffab/r/ancombc.R") or command[1].endswith(
        "microsuite\\diffab\\r\\ancombc.R"
    )
    assert command[-1] == str(tmp_path / "diff.tsv")
    params = captured_params[0]
    assert params["fix_formula"] == "treatment"
    assert params["group"] == "treatment"


@pytest.mark.parametrize("backend", ["aldex2", "maaslin2", "lefse"])
def test_r_diffab_backend_invokes_packaged_script(
    backend: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    commands: list[list[str]] = []

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript")

    def fake_run(command: list[str], **kwargs: object) -> object:
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("subprocess.run", fake_run)

    diff_abundance(
        backend=backend,
        table=table,
        group="treatment",
        output=tmp_path / "diff.tsv",
    )

    assert commands
    command = commands[0]
    assert command[0] == "Rscript"
    assert command[1].endswith(f"microsuite/diffab/r/{backend}.R") or command[1].endswith(
        f"microsuite\\diffab\\r\\{backend}.R"
    )
    assert command[-2:] == ["treatment", str(tmp_path / "diff.tsv")]


def test_diff_abundance_writes_runtime_logs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    run_dir = tmp_path / "run"

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript")
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    diff_abundance(
        backend="ancombc",
        table=table,
        group="treatment",
        output=tmp_path / "diff.tsv",
        run_dir=run_dir,
    )

    assert (run_dir / "command.txt").read_text(encoding="utf-8").startswith("Rscript ")
    assert (run_dir / "stdout.log").read_text(encoding="utf-8") == "ok\n"
    assert (run_dir / "stderr.log").read_text(encoding="utf-8") == ""
    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "diff_abundance"
    assert run["backend"] == "ancombc"
    assert run["inputs"] == {"fix_formula": "treatment", "rand_formula": ""}
    assert run["outputs"] == {"output": str(tmp_path / "diff.tsv")}


@pytest.mark.parametrize("backend", ["aldex2", "maaslin2", "lefse"])
def test_r_diffab_backend_writes_runtime_logs(
    backend: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    run_dir = tmp_path / "run"

    monkeypatch.setattr(shutil, "which", lambda name: "Rscript")
    monkeypatch.setattr(
        "subprocess.run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 0, "ok\n", ""),
    )

    diff_abundance(
        backend=backend,
        table=table,
        group="treatment",
        output=tmp_path / "diff.tsv",
        run_dir=run_dir,
    )

    run = json.loads((run_dir / "run.json").read_text(encoding="utf-8"))
    assert run["task"] == "diff_abundance"
    assert run["backend"] == backend
    assert run["inputs"] == {"group": "treatment"}
    assert run["outputs"] == {"output": str(tmp_path / "diff.tsv")}


@pytest.mark.parametrize("backend", ["aldex2", "maaslin2", "lefse"])
def test_r_diffab_backend_missing_rscript(
    backend: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="Rscript"):
        diff_abundance(
            backend=backend,
            table=table,
            group="treatment",
            output=tmp_path / "diff.tsv",
        )


def test_cli_diff_abundance_help_and_missing_rscript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: None)
    runner = CliRunner()

    result = runner.invoke(
        app,
        [
            "diff_abundance",
            "--backend",
            "ancombc",
            "--table",
            str(table),
            "--group",
            "treatment",
            "-o",
            str(tmp_path / "diff.tsv"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
    assert "Rscript" in str(result.exception)


def test_legacy_diffab_ancombc_command_still_reports_rscript(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    monkeypatch.setattr(shutil, "which", lambda name: None)

    result = CliRunner().invoke(
        app,
        [
            "diffab",
            "ancombc",
            str(table),
            "--group",
            "treatment",
            "-o",
            str(tmp_path / "diff.tsv"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
    assert "Rscript" in str(result.exception)


def test_diff_abundance_threads_runtime_and_image(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = fixture_table(tmp_path)
    captured: dict = {}
    monkeypatch.setattr(
        "microsuite.methods.diff_abundance.run_ancombc",
        lambda adata, **kw: captured.update(kw),
    )
    diff_abundance(
        backend="ancombc",
        table=table,
        group="x",
        output=tmp_path / "o.tsv",
        runtime="docker",
        image="img:2",
    )
    assert captured["runtime"] == "docker" and captured["image"] == "img:2"

    captured_r: dict = {}
    monkeypatch.setattr(
        "microsuite.methods.diff_abundance.run_r_diffab_backend",
        lambda adata, **kw: captured_r.update(kw),
    )
    diff_abundance(
        backend="aldex2",
        table=table,
        group="x",
        output=tmp_path / "o2.tsv",
        runtime="docker",
        image="img:3",
    )
    assert captured_r["runtime"] == "docker" and captured_r["image"] == "img:3"


def test_maaslin2_script_normalizes_library_size() -> None:
    # The caller hands MaAsLin 2 a raw count table. normalization = "NONE" left
    # every result confounded by sequencing depth while still producing a
    # well-formed table of plausible p-values, so nothing surfaced the error.
    # TSS is MaAsLin 2's own default. Shipped wrong in 0.2.0; fixed in 0.2.1.
    script = files("microsuite.diffab.r").joinpath("maaslin2.R").read_text(encoding="utf-8")

    assert 'normalization = "TSS"' in script
    assert 'normalization = "NONE"' not in script


def test_lefse_script_converts_counts_to_relative_abundance() -> None:
    # lefser documents that LEfSe expects relative abundances; handed raw counts
    # it only warns and continues, yielding LDA scores driven by library size.
    # Shipped wrong in 0.2.0; fixed in 0.2.1.
    script = files("microsuite.diffab.r").joinpath("lefse.R").read_text(encoding="utf-8")

    assert "relativeAb" in script
    assert script.index("relativeAb") < script.index("lefser::lefser")
