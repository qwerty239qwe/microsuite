from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.h5ad import read_h5ad, write_h5ad
from microsuite.io.tsv import read_tsv
from microsuite.methods.abundance import abundance_native
from microsuite.methods.normalize import normalize_native
from microsuite.methods.rarefy import rarefy_native
from microsuite.methods.shared_taxa import shared_taxa_native

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def test_normalize_relative_rows_sum_to_one() -> None:
    result = normalize_native(fixture_adata(), method="relative")

    assert np.allclose(np.asarray(result.X).sum(axis=1), 1.0)
    assert result.uns["microsuite_normalize"]["method"] == "relative"


def test_normalize_clr_centers_each_sample() -> None:
    result = normalize_native(fixture_adata(), method="clr", pseudocount=1.0)

    assert np.allclose(np.asarray(result.X).mean(axis=1), 0.0)


def test_prevalence_filter_removes_low_prevalence_features() -> None:
    result = normalize_native(fixture_adata(), method="prevalence-filter", min_prevalence=0.75)

    assert list(result.var_names) == ["f1", "f3", "f4"]


def test_abundance_native_aggregates_taxonomy_level() -> None:
    result = abundance_native(fixture_adata(), level="genus", relative=False)

    assert result.loc["L1S8", "Lactobacillus"] == 10
    assert result.loc["L1S57", "Bacteroides"] == 5


def test_shared_taxa_native_groups_presence() -> None:
    result = shared_taxa_native(fixture_adata(), level="genus", group="body_site")

    lactobacillus = result.set_index("taxon").loc["Lactobacillus"]
    assert lactobacillus["n_groups"] == 2
    assert lactobacillus["groups"] == "gut,tongue"


def test_rarefy_native_fixed_depth_and_seed() -> None:
    result = rarefy_native(fixture_adata(), depth=10, seed=42)

    assert np.asarray(result.X).sum(axis=1).tolist() == [10.0, 10.0, 10.0, 10.0]
    assert result.uns["microsuite_rarefy"]["depth"] == 10


def test_rarefy_native_rejects_low_depth_samples() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="below rarefaction depth"):
        rarefy_native(fixture_adata(), depth=20)


def test_rarefy_native_never_exceeds_original_counts() -> None:
    # Without-replacement subsampling can never draw a feature more times than it
    # was originally observed. A with-replacement (multinomial) draw would.
    adata = fixture_adata()
    original = np.asarray(adata.X, dtype=np.float64)
    for seed in range(50):
        result = np.asarray(rarefy_native(adata, depth=10, seed=seed).X)
        assert np.all(result <= original)


def test_rarefy_native_at_full_depth_returns_table_unchanged() -> None:
    # At depth == sample total, true rarefaction returns the original counts exactly.
    adata = fixture_adata()
    original = np.asarray(adata.X, dtype=np.float64)
    depth = int(original.sum(axis=1).min())
    result = np.asarray(rarefy_native(adata, depth=depth, seed=1).X)
    unchanged = original[original.sum(axis=1) == depth]
    got = result[original.sum(axis=1) == depth]
    assert np.array_equal(got, unchanged)


def test_h5ad_roundtrip_for_normalize(tmp_path: Path) -> None:
    table = tmp_path / "table.h5ad"
    output = tmp_path / "relative.h5ad"
    write_h5ad(fixture_adata(), table)
    write_h5ad(normalize_native(read_h5ad(table), method="relative"), output)

    assert pd.DataFrame(read_h5ad(output).obs).shape[0] == 4


def test_normalize_clr_writes_layer_preserves_x(tmp_path) -> None:
    import numpy as np

    from microsuite.io.h5ad import read_h5ad, write_h5ad
    from microsuite.io.tsv import read_matrix_tsv
    from microsuite.methods.normalize import normalize, normalize_native

    counts = tmp_path / "counts.tsv"
    counts.write_text("feature_id\ts1\ts2\nA\t5\t1\nB\t3\t9\n", encoding="utf-8")

    adata_in = read_matrix_tsv(counts)
    src = tmp_path / "in.h5ad"
    write_h5ad(adata_in, src)

    out = tmp_path / "out.h5ad"
    normalize(backend="native", method="clr", table=src, output=out)
    result = read_h5ad(out)
    # X still raw counts
    assert np.allclose(result.X, adata_in.X)
    # clr stored as a layer, matching normalize_native
    assert "clr" in result.layers
    expected = normalize_native(adata_in, method="clr").X
    assert np.allclose(result.layers["clr"], expected)


def test_normalize_prevalence_filter_filters_features(tmp_path) -> None:
    from microsuite.io.h5ad import read_h5ad, write_h5ad
    from microsuite.io.tsv import read_matrix_tsv
    from microsuite.methods.normalize import normalize

    counts = tmp_path / "counts.tsv"
    # feature B present in only 1 of 2 samples -> prevalence 0.5
    counts.write_text("feature_id\ts1\ts2\nA\t5\t1\nB\t0\t9\n", encoding="utf-8")
    src = tmp_path / "in.h5ad"
    write_h5ad(read_matrix_tsv(counts), src)
    out = tmp_path / "out.h5ad"
    normalize(
        backend="native", method="prevalence-filter", table=src, output=out, min_prevalence=0.75
    )
    result = read_h5ad(out)
    assert list(result.var_names) == ["A"]  # B filtered
    assert "prevalence-filter" not in result.layers  # structural, not a layer
