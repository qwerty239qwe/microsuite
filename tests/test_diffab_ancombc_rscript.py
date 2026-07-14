from __future__ import annotations

import json
import shutil
import subprocess
from importlib.resources import files
from pathlib import Path

import pytest

SCRIPT = str(files("microsuite.diffab.r").joinpath("ancombc.R"))


def _have_r_jsonlite() -> bool:
    if shutil.which("Rscript") is None:
        return False
    r = subprocess.run(
        ["Rscript", "-e", 'quit(status = !requireNamespace("jsonlite", quietly = TRUE))'],
        capture_output=True,
    )
    return r.returncode == 0


pytestmark = pytest.mark.skipif(
    not _have_r_jsonlite(), reason="Rscript with jsonlite not available"
)


def _write_inputs(tmp_path: Path, params: dict) -> tuple[Path, Path, Path, Path]:
    counts = tmp_path / "counts.tsv"
    counts.write_text("\ts1\ts2\ts3\ts4\nF1\t5\t1\t8\t2\nF2\t0\t3\t1\t4\n", encoding="utf-8")
    meta = tmp_path / "meta.tsv"
    # phase and time are perfectly confounded -> rank-deficient when combined
    meta.write_text(
        "\tphase\ttime\tsubject\ns1\tpre\t0\ta\ns2\tpre\t0\tb\ns3\tpost\t7\ta\ns4\tpost\t7\tb\n",
        encoding="utf-8",
    )
    pj = tmp_path / "params.json"
    pj.write_text(json.dumps(params), encoding="utf-8")
    out = tmp_path / "out.tsv"
    return counts, meta, pj, out


def test_ancombc_r_rejects_rank_deficient_design(tmp_path: Path) -> None:
    counts, meta, pj, out = _write_inputs(
        tmp_path, {"fix_formula": "phase + time", "reference": {}}
    )
    r = subprocess.run(
        ["Rscript", SCRIPT, str(counts), str(meta), str(pj), str(out)],
        capture_output=True,
        text=True,
    )
    assert r.returncode != 0
    assert "rank" in (r.stdout + r.stderr).lower()


def test_ancombc_r_reports_missing_formula_column(tmp_path: Path) -> None:
    counts, meta, pj, out = _write_inputs(tmp_path, {"fix_formula": "nope", "reference": {}})
    r = subprocess.run(
        ["Rscript", SCRIPT, str(counts), str(meta), str(pj), str(out)],
        capture_output=True,
        text=True,
    )
    assert r.returncode != 0
    assert "nope" in (r.stdout + r.stderr)
