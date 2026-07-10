from __future__ import annotations

import json
from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.dada2_manifest import (
    MANIFEST_FILENAME,
    build_manifest,
    read_r_params,
    write_manifest,
)

FIXTURE = Path(__file__).parent / "fixtures" / "dada2_r_params_paired.json"


def _wrapper() -> dict:
    return {
        "microsuite_version": "9.9.9",
        "backend": "dada2-r",
        "runtime": "docker",
        "image": "ghcr.io/example/r-dada2:latest",
        "mode": "paired",
        "paired": True,
        "threads": 4,
        "input_dir": "/data/reads",
        "output_table": "/out/table.tsv",
        "output_rep_seqs": "/out/rep.fasta",
        "output_stats": "/out/stats.tsv",
        "output_plot_dir": None,
        "created_at": "2026-07-10T00:00:00+00:00",
        "command": "Rscript dada2_denoise.R --paired",
    }


def test_read_r_params_ok() -> None:
    params = read_r_params(FIXTURE)
    assert params["min_overlap"] == 12
    assert params["dada2_version"] == "1.30.0"
    assert params["trim_left"] is None  # single-only key, null for paired run


def test_read_r_params_missing(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        read_r_params(tmp_path / "nope.json")


def test_read_r_params_malformed(tmp_path: Path) -> None:
    bad = tmp_path / "bad.json"
    bad.write_text("{not json", encoding="utf-8")
    with pytest.raises(MicrobiomeSuiteError):
        read_r_params(bad)


def test_build_manifest_splits_and_drops_nulls() -> None:
    manifest = build_manifest(read_r_params(FIXTURE), _wrapper())
    assert manifest["tool"] == {
        "dada2_version": "1.30.0",
        "r_version": "R version 4.3.2 (2023-10-31)",
    }
    dp = manifest["dada2_params"]
    assert dp["min_overlap"] == 12  # resolved default present, not absent
    assert dp["mode"] == "paired"
    assert "trim_left" not in dp  # single-only null key dropped
    assert "dada2_version" not in dp  # versions live under tool, not params
    assert manifest["run"]["runtime"] == "docker"
    assert manifest["run"]["image"] == "ghcr.io/example/r-dada2:latest"


def test_write_manifest_roundtrip(tmp_path: Path) -> None:
    manifest = build_manifest(read_r_params(FIXTURE), _wrapper())
    path = write_manifest(manifest, tmp_path)
    assert path.name == MANIFEST_FILENAME
    assert json.loads(path.read_text())["dada2_params"]["min_overlap"] == 12
