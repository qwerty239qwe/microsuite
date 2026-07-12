from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.methods.dada2_sweep import GridPoint, run_dada2_sweep

STATS = (
    "\tinput\tfiltered\tdenoised_f\tdenoised_r\tmerged\tnonchim\n"
    "s1\t1000\t900\t880\t870\t800\t700\n"
)


def _fake_denoise_factory(fail_names=()):
    def fake_denoise(*, backend, demux, output_table, output_rep_seqs, output_stats, **kw):
        name = Path(output_table).parent.name
        if name in fail_names:
            raise MicrobiomeSuiteError(f"denoise failed for {name}")
        Path(output_table).write_text("\ts1\nASV1\t5\n", encoding="utf-8")
        Path(output_rep_seqs).write_text(">ASV1\nAAA\n", encoding="utf-8")
        Path(output_stats).write_text(STATS, encoding="utf-8")

    return fake_denoise


def test_run_dada2_sweep_writes_summary(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr("microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory())
    grid = [
        GridPoint(name="baseline", params={"max_ee_f": 2}, is_baseline=True),
        GridPoint(name="relaxed", params={"max_ee_f": 3}, is_baseline=False),
    ]
    out = run_dada2_sweep(
        input_dir=tmp_path / "reads", mode="paired", output_dir=tmp_path / "out", grid=grid
    )
    assert out.exists()
    df = pd.read_csv(out, sep="\t")
    assert list(df["name"]) == ["baseline", "relaxed"]
    assert set(df["status"]) == {"ok"}


def test_run_dada2_sweep_failed_point_recorded(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        "microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory(fail_names={"relaxed"})
    )
    grid = [
        GridPoint(name="baseline", params={"max_ee_f": 2}, is_baseline=True),
        GridPoint(name="relaxed", params={"max_ee_f": 1}, is_baseline=False),
    ]
    with pytest.warns(UserWarning, match="relaxed"):
        out = run_dada2_sweep(
            input_dir=tmp_path / "reads", mode="paired", output_dir=tmp_path / "out", grid=grid
        )
    df = pd.read_csv(out, sep="\t")
    assert df[df["name"] == "relaxed"].iloc[0]["status"] == "failed"


def test_run_dada2_sweep_baseline_failure_raises(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        "microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory(fail_names={"baseline"})
    )
    grid = [GridPoint(name="baseline", params={}, is_baseline=True)]
    with pytest.raises(MicrobiomeSuiteError, match="[Bb]aseline"):
        run_dada2_sweep(
            input_dir=tmp_path / "reads", mode="paired", output_dir=tmp_path / "out", grid=grid
        )


def test_cli_denoise_sweep_axes(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr("microsuite.methods.dada2_sweep.denoise", _fake_denoise_factory())
    (tmp_path / "reads").mkdir()
    r = CliRunner().invoke(
        app,
        [
            "denoise-sweep",
            "--input-dir",
            str(tmp_path / "reads"),
            "--mode",
            "paired",
            "--output-dir",
            str(tmp_path / "out"),
            "--max-ee-f",
            "2,3",
        ],
    )
    assert r.exit_code == 0, r.stdout
    assert (tmp_path / "out" / "dada2_sweep_summary.tsv").exists()


def test_cli_denoise_sweep_both_sources_errors(tmp_path) -> None:
    cfg = tmp_path / "g.json"
    cfg.write_text("[]", encoding="utf-8")
    r = CliRunner().invoke(
        app,
        [
            "denoise-sweep",
            "--input-dir",
            str(tmp_path / "reads"),
            "--mode",
            "paired",
            "--output-dir",
            str(tmp_path / "out"),
            "--grid-config",
            str(cfg),
            "--max-ee-f",
            "2",
        ],
    )
    assert r.exit_code != 0
