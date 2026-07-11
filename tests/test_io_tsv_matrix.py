from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.tsv import read_count_matrix, read_matrix_tsv


def _write(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_read_count_matrix_basic(tmp_path: Path) -> None:
    p = _write(tmp_path / "c.tsv", "feature_id\ts1\ts2\nASV1\t5\t1\nASV2\t0\t3\n")
    m = read_count_matrix(p)
    assert list(m.index) == ["ASV1", "ASV2"]
    assert list(m.columns) == ["s1", "s2"]
    assert m.index.name == "feature_id"
    assert m.loc["ASV1", "s2"] == 1


def test_read_count_matrix_sanitizes_rank_first_column(tmp_path: Path) -> None:
    p = _write(tmp_path / "g.tsv", "genus\ts1\ts2\nBacteroides\t5\t1\nPrevotella\t0\t3\n")
    with pytest.warns(UserWarning, match="feature_id"):
        m = read_count_matrix(p)
    assert m.index.name == "feature_id"
    assert list(m.index) == ["Bacteroides", "Prevotella"]  # IDs preserved


def test_read_count_matrix_no_warn_normal_header(tmp_path: Path, recwarn) -> None:
    p = _write(tmp_path / "c.tsv", "#OTU ID\ts1\nASV1\t5\n")
    read_count_matrix(p)
    assert not any(issubclass(w.category, UserWarning) for w in recwarn.list)


def test_read_count_matrix_rejects_empty_and_dups(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        read_count_matrix(_write(tmp_path / "e.tsv", "feature_id\n"))
    with pytest.raises(MicrobiomeSuiteError):
        read_count_matrix(_write(tmp_path / "d.tsv", "feature_id\ts1\nASV1\t5\nASV1\t3\n"))


def test_read_matrix_tsv_no_metadata(tmp_path: Path) -> None:
    p = _write(tmp_path / "c.tsv", "feature_id\ts1\ts2\nASV1\t5\t1\nASV2\t0\t3\n")
    adata = read_matrix_tsv(p)
    assert adata.shape == (2, 2)  # 2 samples (obs) x 2 features (var)
    assert list(adata.obs_names) == ["s1", "s2"]
    assert list(adata.var_names) == ["ASV1", "ASV2"]
    assert adata.obs.shape[1] == 0  # empty metadata


def test_read_tsv_with_rank_named_first_column(tmp_path: Path) -> None:
    from microsuite.io.h5ad import write_h5ad
    from microsuite.io.tsv import read_tsv

    _write(tmp_path / "g.tsv", "genus\ts1\ts2\nBacteroides\t5\t1\nPrevotella\t0\t3\n")
    _write(tmp_path / "m.tsv", "sample\tgroup\ns1\ta\ns2\tb\n")
    with pytest.warns(UserWarning, match="feature_id"):
        adata = read_tsv(tmp_path / "g.tsv", tmp_path / "m.tsv")
    # the previously-breaking write must now succeed
    write_h5ad(adata, tmp_path / "out.h5ad")
    assert (tmp_path / "out.h5ad").exists()
    assert list(adata.var_names) == ["Bacteroides", "Prevotella"]
