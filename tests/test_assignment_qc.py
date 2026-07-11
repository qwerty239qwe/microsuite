from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.taxonomy import LEVELS
from microsuite.methods.assignment_qc import (
    POOLED_LABEL,
    deepest_rank_distribution,
    summarize_assignment,
    write_assignment_summary,
)


def _fixture() -> ad.AnnData:
    # F1 assigned to species (all ranks), F2 assigned to genus only, F3 unassigned.
    rank_values = {
        "kingdom": ["Bacteria", "Bacteria", ""],
        "phylum": ["Firmicutes", "Bacteroidetes", ""],
        "class": ["Bacilli", "Bacteroidia", ""],
        "order": ["Lactobacillales", "Bacteroidales", ""],
        "family": ["Lactobacillaceae", "Prevotellaceae", ""],
        "genus": ["Lactobacillus", "Prevotella", ""],
        "species": ["L. casei", "", ""],
    }
    var = pd.DataFrame(rank_values, index=["F1", "F2", "F3"])
    # samples s1 (all present), s2 (F1 absent)
    X = np.array([[10.0, 5.0, 1.0], [0.0, 3.0, 2.0]])
    return ad.AnnData(X=X, obs=pd.DataFrame(index=["s1", "s2"]), var=var)


def test_summarize_assignment_overall_counts() -> None:
    df = summarize_assignment(_fixture())
    pooled = df[df["sample"] == POOLED_LABEL].set_index("rank")
    # species: only F1 assigned -> 1 assigned, 2 unassigned features
    assert pooled.loc["species", "assigned_features"] == 1
    assert pooled.loc["species", "unassigned_features"] == 2
    # genus: F1 + F2 assigned -> 2 assigned, 1 unassigned
    assert pooled.loc["genus", "assigned_features"] == 2
    assert pooled.loc["genus", "unassigned_features"] == 1
    # pooled reads at species: assigned = F1 reads (10) ; unassigned = F2+F3 (5+3+1+2)=11
    assert pooled.loc["species", "assigned_reads"] == 10.0
    assert pooled.loc["species", "unassigned_reads"] == 11.0
    # assigned_features is non-increasing kingdom -> species (nested assignment)
    ranks = [r for r in LEVELS]
    seq = [pooled.loc[r, "assigned_features"] for r in ranks]
    assert all(seq[i] >= seq[i + 1] for i in range(len(seq) - 1))


def test_summarize_assignment_per_sample_presence() -> None:
    df = summarize_assignment(_fixture())
    s2 = df[df["sample"] == "s2"].set_index("rank")
    # in s2, F1 absent (0 reads); present = {F2, F3}
    # genus: F2 assigned, F3 not -> 1 assigned, 1 unassigned
    assert s2.loc["genus", "assigned_features"] == 1
    assert s2.loc["genus", "unassigned_features"] == 1
    # s2 reads: total 5? no -> F2=3, F3=2 => genus assigned_reads=3, unassigned=2
    assert s2.loc["genus", "assigned_reads"] == 3.0
    assert s2.loc["genus", "unassigned_reads"] == 2.0
    assert s2.loc["genus", "assigned_read_frac"] == pytest.approx(3.0 / 5.0)


def test_deepest_rank_distribution() -> None:
    dist = deepest_rank_distribution(_fixture())
    assert dist["species"] == 1  # F1
    assert dist["genus"] == 1  # F2
    assert dist["Unassigned"] == 1  # F3
    assert int(dist.sum()) == 3  # each ASV counted once


def _fixture_with_nan() -> ad.AnnData:
    # F1 assigned to species (all ranks), F2 assigned to genus only (with a
    # real NaN at species, simulating an uncovered feature from
    # taxonomy.reindex), F3 unassigned via NaN at every rank.
    rank_values = {
        "kingdom": ["Bacteria", "Bacteria", np.nan],
        "phylum": ["Firmicutes", "Bacteroidetes", np.nan],
        "class": ["Bacilli", "Bacteroidia", np.nan],
        "order": ["Lactobacillales", "Bacteroidales", np.nan],
        "family": ["Lactobacillaceae", "Prevotellaceae", np.nan],
        "genus": ["Lactobacillus", "Prevotella", np.nan],
        "species": ["L. casei", np.nan, np.nan],
    }
    var = pd.DataFrame(rank_values, index=["F1", "F2", "F3"])
    X = np.array([[10.0, 5.0, 1.0], [0.0, 3.0, 2.0]])
    return ad.AnnData(X=X, obs=pd.DataFrame(index=["s1", "s2"]), var=var)


def test_summarize_assignment_treats_nan_as_unassigned() -> None:
    # Regression test: a real np.nan (e.g. from taxonomy.reindex on features
    # missing from the taxonomy table) must count as UNASSIGNED, not
    # assigned. Previously `.astype(str)` turned NaN into the string "nan",
    # which is != "" and so was incorrectly counted as assigned.
    df = summarize_assignment(_fixture_with_nan())
    pooled = df[df["sample"] == POOLED_LABEL].set_index("rank")
    # species: only F1 assigned (F2, F3 have NaN) -> 1 assigned, 2 unassigned
    assert pooled.loc["species", "assigned_features"] == 1
    assert pooled.loc["species", "unassigned_features"] == 2
    # genus: F1 + F2 assigned, F3 has NaN -> 2 assigned, 1 unassigned
    assert pooled.loc["genus", "assigned_features"] == 2
    assert pooled.loc["genus", "unassigned_features"] == 1


def test_deepest_rank_distribution_treats_nan_as_unassigned() -> None:
    dist = deepest_rank_distribution(_fixture_with_nan())
    assert dist["species"] == 1  # F1
    assert dist["genus"] == 1  # F2 (NaN tail from species onward)
    assert dist["Unassigned"] == 1  # F3 (all-NaN)
    assert int(dist.sum()) == 3


def test_no_taxonomy_ranks_raises() -> None:
    adata = ad.AnnData(
        X=np.array([[1.0, 2.0]]),
        obs=pd.DataFrame(index=["s1"]),
        var=pd.DataFrame(index=["F1", "F2"]),
    )
    with pytest.raises(MicrobiomeSuiteError):
        summarize_assignment(adata)


def test_write_assignment_summary_roundtrip(tmp_path: Path) -> None:
    df = summarize_assignment(_fixture())
    out = write_assignment_summary(df, tmp_path / "s.tsv")
    assert out.exists()
    back = pd.read_csv(out, sep="\t")
    assert set(["sample", "rank", "assigned_features", "assigned_read_frac"]).issubset(back.columns)
    assert (back["sample"] == POOLED_LABEL).any()
