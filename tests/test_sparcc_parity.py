from __future__ import annotations

from functools import lru_cache
from pathlib import Path

import numpy as np
import pytest

from microsuite.methods._sparcc import SparCCResult, estimate_sparcc

FIXTURE_DIR = Path(__file__).parent / "fixtures" / "sparcc"
REFERENCE_SEEDS = (10010, 10011, 10012)
ALLOWED_MAE = 0.05300942276035425
OLD_CLR_REFERENCE_MAE = {
    "dense": 0.09383232047077433,
    "zero": 0.09221407177885772,
}


def _read_matrix(filename: str) -> np.ndarray:
    return np.loadtxt(
        FIXTURE_DIR / filename,
        delimiter="\t",
        skiprows=1,
        usecols=range(1, 11),
    )


def _off_diagonal_mae(left: np.ndarray, right: np.ndarray) -> float:
    upper = np.triu_indices_from(left, k=1)
    return float(np.mean(np.abs(left[upper] - right[upper])))


def _clr_pearson(counts: np.ndarray) -> np.ndarray:
    logged = np.log(counts + 1.0)
    clr_values = logged - logged.mean(axis=1, keepdims=True)
    return np.corrcoef(clr_values, rowvar=False)


@lru_cache(maxsize=2)
def _case(dataset: str) -> tuple[np.ndarray, SparCCResult, np.ndarray, np.ndarray]:
    counts = _read_matrix(f"{dataset}_counts.tsv")
    result = estimate_sparcc(counts, seed=10010)
    references = np.stack(
        [_read_matrix(f"{dataset}_reference_cor_seed_{seed}.tsv") for seed in REFERENCE_SEEDS]
    )
    reference_median = np.median(references, axis=0)
    truth = _read_matrix("latent_correlation.tsv")
    return counts, result, reference_median, truth


@pytest.mark.parametrize("dataset", ["dense", "zero"])
def test_outer_estimator_matches_frozen_reference_variability(dataset: str) -> None:
    _, result, reference, _ = _case(dataset)

    assert _off_diagonal_mae(result.correlation, reference) <= ALLOWED_MAE


@pytest.mark.parametrize("dataset", ["dense", "zero"])
def test_outer_estimator_is_materially_closer_than_old_clr_backend(dataset: str) -> None:
    _, result, reference, _ = _case(dataset)

    assert _off_diagonal_mae(result.correlation, reference) <= (
        OLD_CLR_REFERENCE_MAE[dataset] / 3.0
    )


@pytest.mark.parametrize("dataset", ["dense", "zero"])
def test_outer_estimator_improves_truth_mae_and_recovers_strong_signs(dataset: str) -> None:
    counts, result, _, truth = _case(dataset)
    clr_correlation = _clr_pearson(counts)
    sparcc_truth_mae = _off_diagonal_mae(result.correlation, truth)
    clr_truth_mae = _off_diagonal_mae(clr_correlation, truth)

    assert sparcc_truth_mae <= 0.9 * clr_truth_mae
    strong_pairs = np.argwhere(np.triu(np.abs(truth) >= 0.7, k=1))
    assert strong_pairs.shape == (3, 2)
    for left, right in strong_pairs:
        assert np.sign(result.correlation[left, right]) == np.sign(truth[left, right])
