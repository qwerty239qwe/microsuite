from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_count_matrix, read_matrix_tsv
from microsuite.methods.normalize import normalize_native
from microsuite.methods.table_io import export_table, normalize_table


def _counts(tmp_path: Path) -> Path:
    p = tmp_path / "counts.tsv"
    p.write_text("feature_id\ts1\ts2\nA\t5\t1\nB\t3\t9\n", encoding="utf-8")
    return p


def test_export_table_roundtrip(tmp_path: Path) -> None:
    src = tmp_path / "in.h5ad"
    write_h5ad(read_matrix_tsv(_counts(tmp_path)), src)
    out = tmp_path / "x.tsv"
    export_table(table=src, output=out)
    m = read_count_matrix(out)
    assert list(m.index) == ["A", "B"]
    assert list(m.columns) == ["s1", "s2"]
    assert m.loc["B", "s2"] == 9


def test_export_table_layer_and_metadata(tmp_path: Path) -> None:
    adata = read_matrix_tsv(_counts(tmp_path))
    adata.layers["clr"] = normalize_native(adata, method="clr").X
    adata.obs["group"] = ["a", "b"]
    src = tmp_path / "in.h5ad"
    write_h5ad(adata, src)
    out = tmp_path / "clr.tsv"
    meta = tmp_path / "meta.tsv"
    export_table(table=src, output=out, layer="clr", metadata=meta)
    exported = read_count_matrix(out)
    assert np.allclose(exported.to_numpy(), adata.layers["clr"].T)
    assert meta.exists() and "group" in meta.read_text()


def test_export_table_missing_layer_errors(tmp_path: Path) -> None:
    src = tmp_path / "in.h5ad"
    write_h5ad(read_matrix_tsv(_counts(tmp_path)), src)
    with pytest.raises(MicrobiomeSuiteError, match="layer"):
        export_table(table=src, output=tmp_path / "x.tsv", layer="nope")


def test_normalize_table_clr_tsv_to_tsv(tmp_path: Path) -> None:
    out = tmp_path / "clr.tsv"
    normalize_table(method="clr", input_path=_counts(tmp_path), output=out)
    exported = read_count_matrix(out)
    expected = normalize_native(read_matrix_tsv(_counts(tmp_path)), method="clr").X
    assert np.allclose(exported.to_numpy(), expected.T)


def test_table_cli_smoke(tmp_path: Path) -> None:
    runner = CliRunner()
    out = tmp_path / "clr.tsv"
    result = runner.invoke(
        app,
        [
            "table",
            "normalize",
            "--method",
            "clr",
            "--input",
            str(_counts(tmp_path)),
            "-o",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert out.exists()
