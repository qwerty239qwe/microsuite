from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from microsuite._errors import MicrobiomeSuiteError


@dataclass(frozen=True)
class SparCCResult:
    covariance: np.ndarray
    correlation: np.ndarray


def estimate_sparcc(
    counts: np.ndarray,
    *,
    iterations: int = 20,
    inner_iterations: int = 10,
    exclusion_threshold: float = 0.1,
    pseudocount: float = 1.0,
    seed: int = 0,
) -> SparCCResult:
    """Estimate SparCC covariance and correlation matrices from raw counts."""
    validated_counts = _validate_inputs(
        counts,
        iterations=iterations,
        inner_iterations=inner_iterations,
        exclusion_threshold=exclusion_threshold,
        pseudocount=pseudocount,
        seed=seed,
    )
    rng = np.random.default_rng(seed)
    _dirichlet_normalize(validated_counts, pseudocount=pseudocount, rng=rng)
    raise NotImplementedError("SparCC basis algebra is implemented in the next plan task.")


def _validate_inputs(
    counts: np.ndarray,
    *,
    iterations: int,
    inner_iterations: int,
    exclusion_threshold: float,
    pseudocount: float,
    seed: int,
) -> np.ndarray:
    validated_counts = _validate_counts(counts)
    _validate_positive_integer(iterations, name="iterations")
    _validate_positive_integer(inner_iterations, name="inner_iterations")
    _validate_threshold(exclusion_threshold)
    _validate_pseudocount(pseudocount)
    _validate_seed(seed)
    return validated_counts


def _validate_counts(counts: np.ndarray) -> np.ndarray:
    try:
        array = np.asarray(counts)
    except (TypeError, ValueError) as exc:
        raise MicrobiomeSuiteError(
            "SparCC counts must be a numeric two-dimensional array."
        ) from exc

    if array.ndim != 2:
        raise MicrobiomeSuiteError("SparCC counts must be a two-dimensional array.")
    if array.shape[0] < 2:
        raise MicrobiomeSuiteError("SparCC requires at least two samples.")
    if array.shape[1] < 3:
        raise MicrobiomeSuiteError("SparCC requires at least three features.")
    if array.dtype.kind not in "iuf":
        raise MicrobiomeSuiteError("SparCC counts must contain only real numeric values.")

    try:
        numeric = array.astype(np.float64, copy=True)
    except (TypeError, ValueError, OverflowError) as exc:
        raise MicrobiomeSuiteError("SparCC counts must contain only numeric values.") from exc

    if not np.isfinite(numeric).all():
        raise MicrobiomeSuiteError("SparCC counts must contain only finite values.")
    if (numeric < 0).any():
        raise MicrobiomeSuiteError("SparCC counts must be nonnegative.")
    if not np.equal(numeric, np.floor(numeric)).all():
        raise MicrobiomeSuiteError("SparCC counts must be integer-valued.")
    if (numeric.sum(axis=1) == 0).any():
        raise MicrobiomeSuiteError("SparCC counts contain an all-zero sample.")
    if (numeric.sum(axis=0) == 0).any():
        raise MicrobiomeSuiteError("SparCC counts contain an all-zero feature.")
    return numeric


def _validate_positive_integer(value: int, *, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 1:
        raise MicrobiomeSuiteError(f"SparCC {name} must be an integer greater than zero.")


def _validate_threshold(value: float) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not np.isfinite(value)
        or not 0.0 <= value <= 1.0
    ):
        raise MicrobiomeSuiteError("SparCC exclusion_threshold must be between zero and one.")


def _validate_pseudocount(value: float) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, Real)
        or not np.isfinite(value)
        or value <= 0.0
    ):
        raise MicrobiomeSuiteError("SparCC pseudocount must be finite and greater than zero.")


def _validate_seed(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, Integral) or value < 0:
        raise MicrobiomeSuiteError("SparCC seed must be a nonnegative integer.")


def _dirichlet_normalize(
    counts: np.ndarray,
    *,
    pseudocount: float,
    rng: np.random.Generator,
) -> np.ndarray:
    """Draw one Dirichlet composition per sample using vectorized gamma draws."""
    concentrations = counts + pseudocount
    draws = rng.gamma(shape=concentrations, scale=1.0)
    totals = draws.sum(axis=1, keepdims=True)
    if not np.isfinite(draws).all() or (draws <= 0.0).any() or not np.isfinite(totals).all():
        raise MicrobiomeSuiteError("SparCC Dirichlet normalization produced invalid values.")
    return draws / totals
