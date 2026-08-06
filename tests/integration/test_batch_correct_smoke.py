"""Real-execution smoke test for the batch-correction backends.

Gated behind ``MICROSUITE_RUN_BATCH_SMOKE=1`` and a container engine, because
it runs the real R backends and the vegan image.

**Why this test exists.** Every other batch test mocks the subprocess, so
together they prove only that we build the commands we meant to build. The
mothur work shipped fourteen defects past a suite of exactly that kind, each
producing a complete, well-formed, wrong result.

The assertion is deliberately a *pair*. A backend that shrinks the batch effect
by flattening every difference in the table scores perfectly on the batch term
alone. Both terms are checked, always.

The dataset is generated deterministically rather than committed, so there is
no fixture to drift.
"""

from __future__ import annotations

import os
import shutil

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite.batch.correct import run_batch_correction
from microsuite.diversity.beta import beta_diversity
from microsuite.diversity.ecology import beta_significance

pytestmark = [
    pytest.mark.skipif(
        os.environ.get("MICROSUITE_RUN_BATCH_SMOKE") != "1",
        reason="set MICROSUITE_RUN_BATCH_SMOKE=1 to run the real batch backends",
    ),
    pytest.mark.skipif(shutil.which("docker") is None, reason="docker is not installed"),
]

N_PER_CELL = 8
N_FEATURES = 60


def _two_batch_dataset(seed: int = 0) -> ad.AnnData:
    """Two batches crossed with two groups, each carrying a distinct shift.

    The design is balanced: every batch contains both groups, so batch and group
    are separable. A correction that removes the group effect along with the
    batch effect is a failure, and this layout is what makes that visible.
    """
    rng = np.random.default_rng(seed)
    rows, batches, groups, names = [], [], [], []
    batch_shift = np.ones(N_FEATURES)
    batch_shift[: N_FEATURES // 2] = 4.0  # multiplicative, half the features
    group_shift = np.ones(N_FEATURES)
    group_shift[N_FEATURES // 2 :] = 2.0  # a smaller effect on the other half

    for batch in ("A", "B"):
        for group in ("case", "ctrl"):
            for replicate in range(N_PER_CELL):
                base = rng.gamma(shape=2.0, scale=30.0, size=N_FEATURES)
                if batch == "B":
                    base = base * batch_shift
                if group == "case":
                    base = base * group_shift
                rows.append(rng.poisson(base).astype(float))
                batches.append(batch)
                groups.append(group)
                names.append(f"{batch}{group}{replicate}")

    adata = ad.AnnData(np.vstack(rows))
    adata.obs_names = names
    adata.var_names = [f"f{i}" for i in range(N_FEATURES)]
    adata.obs = pd.DataFrame({"run_id": batches, "group": groups}, index=names)
    return adata


def _batch_and_group_r2(adata: ad.AnnData) -> tuple[float, float]:
    """PERMANOVA variance explained by batch and by group, via vegan adonis2."""
    distances = beta_diversity(adata, "bray-curtis")
    result = beta_significance(
        distances,
        pd.DataFrame(adata.obs),
        method="adonis2",
        formula="run_id + group",
        backend="vegan",
        permutations=199,
        seed=0,
        runtime="docker",
    )
    indexed = result.set_index("term")["r_squared"]
    return float(indexed["run_id"]), float(indexed["group"])


def _batch_and_group_r2_euclidean(adata: ad.AnnData, *, clr: bool = False) -> tuple[float, float]:
    """Same partition, on Euclidean distance -- the right geometry for CLR data.

    Bray-Curtis on CLR values is meaningless (the values are signed), so
    comparing a CLR-space correction against a Bray-Curtis baseline would
    compare two different quantities and call the difference an improvement.
    """
    from scipy.spatial.distance import pdist, squareform

    from microsuite.methods.normalize import normalize_native

    source = normalize_native(adata, method="clr") if clr else adata
    values = np.asarray(source.X, dtype=float)
    distances = pd.DataFrame(
        squareform(pdist(values, metric="euclidean")),
        index=list(map(str, adata.obs_names)),
        columns=list(map(str, adata.obs_names)),
    )
    result = beta_significance(
        distances,
        pd.DataFrame(adata.obs),
        method="adonis2",
        formula="run_id + group",
        backend="vegan",
        permutations=199,
        seed=0,
        runtime="docker",
    )
    indexed = result.set_index("term")["r_squared"]
    return float(indexed["run_id"]), float(indexed["group"])


def test_uncorrected_dataset_has_the_effects_the_test_assumes() -> None:
    # If the generator stops producing a batch effect, every downstream
    # assertion below becomes vacuously true. Check the premise first.
    batch_r2, group_r2 = _batch_and_group_r2(_two_batch_dataset())
    assert batch_r2 > 0.10, f"generated batch effect too small: {batch_r2}"
    assert group_r2 > 0.02, f"generated group effect too small: {group_r2}"


def test_mmuphin_shrinks_batch_and_keeps_group() -> None:
    adata = _two_batch_dataset()
    before_batch, before_group = _batch_and_group_r2(adata)

    corrected = run_batch_correction(
        adata, backend="mmuphin", batch="run_id", covariates=["group"], runtime="docker"
    )
    after_batch, after_group = _batch_and_group_r2(corrected)

    assert after_batch < before_batch * 0.5, (
        f"batch R2 did not shrink: {before_batch:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group * 0.5, (
        f"biological signal was flattened along with the batch effect: "
        f"{before_group:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "relative"


def test_combat_seq_shrinks_batch_keeps_group_and_returns_counts() -> None:
    adata = _two_batch_dataset()
    before_batch, before_group = _batch_and_group_r2(adata)

    corrected = run_batch_correction(
        adata, backend="combat-seq", batch="run_id", covariates=["group"], runtime="docker"
    )
    after_batch, after_group = _batch_and_group_r2(corrected)

    assert after_batch < before_batch * 0.5, (
        f"batch R2 did not shrink: {before_batch:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group * 0.5, (
        f"biological signal was flattened: {before_group:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "counts"
    # The whole point of this backend is that ANCOM-BC can consume its output.
    values = np.asarray(corrected.X)
    assert np.allclose(values, np.round(values)), "ComBat_seq returned non-integer counts"


def test_conqur_shrinks_batch_and_keeps_group() -> None:
    adata = _two_batch_dataset()
    before_batch, before_group = _batch_and_group_r2(adata)

    corrected = run_batch_correction(
        adata, backend="conqur", batch="run_id", covariates=["group"], runtime="docker"
    )
    after_batch, after_group = _batch_and_group_r2(corrected)

    assert after_batch < before_batch * 0.5, (
        f"batch R2 did not shrink: {before_batch:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group * 0.5, (
        f"biological signal was flattened: {before_group:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "counts"
    values = np.asarray(corrected.X)
    assert np.allclose(values, np.round(values)), "ConQuR returned non-integer counts"


def test_plsda_batch_shrinks_batch_keeps_group_and_returns_clr() -> None:
    adata = _two_batch_dataset()

    corrected = run_batch_correction(
        adata, backend="plsda-batch", batch="run_id", target="group", runtime="docker"
    )
    # CLR output is not a distance-compatible abundance table for Bray-Curtis, so
    # this backend is assessed on Euclidean distance over the CLR values, which is
    # the Aitchison distance the method itself is defined against.
    after_batch, after_group = _batch_and_group_r2_euclidean(corrected)
    before_batch_e, before_group_e = _batch_and_group_r2_euclidean(adata, clr=True)

    assert after_batch < before_batch_e * 0.5, (
        f"batch R2 did not shrink: {before_batch_e:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group_e * 0.5, (
        f"biological signal was flattened: {before_group_e:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "clr"


def test_metadict_shrinks_batch_and_keeps_group() -> None:
    adata = _two_batch_dataset()
    before_batch, before_group = _batch_and_group_r2(adata)

    corrected = run_batch_correction(
        adata, backend="metadict", batch="run_id", covariates=["group"], runtime="docker"
    )
    after_batch, after_group = _batch_and_group_r2(corrected)

    assert after_batch < before_batch * 0.5, (
        f"batch R2 did not shrink: {before_batch:.3f} -> {after_batch:.3f}"
    )
    assert after_group > before_group * 0.5, (
        f"biological signal was flattened along with the batch effect: "
        f"{before_group:.3f} -> {after_group:.3f}"
    )
    assert corrected.uns["microsuite"]["value_type"] == "relative"
