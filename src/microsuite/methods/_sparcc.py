from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral, Real

import numpy as np

from microsuite._errors import MicrobiomeSuiteError

_MIN_BASIS_VARIANCE = 1e-4
_CORRELATION_ROUNDOFF_TOLERANCE = 1e-12


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
    raise NotImplementedError(
        "SparCC iterative exclusion and outer aggregation are not implemented yet."
    )


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


def _clr(compositions: np.ndarray) -> np.ndarray:
    """Return sample-wise centered log-ratios for positive compositions."""
    if compositions.ndim != 2 or not np.isfinite(compositions).all():
        raise MicrobiomeSuiteError("SparCC compositions must be a finite two-dimensional array.")
    if (compositions <= 0.0).any():
        raise MicrobiomeSuiteError("SparCC compositions must be strictly positive.")
    logged = np.log(compositions)
    return logged - logged.mean(axis=1, keepdims=True)


def _variation_matrix(clr_values: np.ndarray) -> np.ndarray:
    """Build the Aitchison variation matrix using sample covariance (ddof=1)."""
    covariance = np.cov(clr_values, rowvar=False, ddof=1)
    variances = np.diag(covariance)
    variation = variances[:, None] + variances[None, :] - 2.0 * covariance
    variation = (variation + variation.T) / 2.0
    np.fill_diagonal(variation, 0.0)
    return variation


def _basis_system(variation: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Construct the initial sparse-correlation basis-variance system."""
    feature_count = variation.shape[0]
    coefficients = np.ones((feature_count, feature_count), dtype=np.float64)
    coefficients += np.eye(feature_count) * (feature_count - 2)
    return coefficients, variation.sum(axis=1)


def _solve_basis_variances(coefficients: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Solve basis variances, using a generalized inverse for singular systems."""
    try:
        variances = np.linalg.solve(coefficients, target)
    except np.linalg.LinAlgError:
        try:
            variances = np.linalg.pinv(coefficients) @ target
        except np.linalg.LinAlgError as exc:
            raise MicrobiomeSuiteError("SparCC could not solve the basis system.") from exc
    if not np.isfinite(variances).all():
        raise MicrobiomeSuiteError("SparCC basis variances are not finite.")
    return np.maximum(variances, _MIN_BASIS_VARIANCE)


def _reconstruct_covariance(variation: np.ndarray, variances: np.ndarray) -> np.ndarray:
    """Reconstruct and symmetrize covariance from variation and basis variance."""
    covariance = 0.5 * (variances[:, None] + variances[None, :] - variation)
    covariance = (covariance + covariance.T) / 2.0
    np.fill_diagonal(covariance, variances)
    if not np.isfinite(covariance).all():
        raise MicrobiomeSuiteError("SparCC reconstructed covariance is not finite.")
    return covariance


def _covariance_to_correlation(covariance: np.ndarray) -> np.ndarray:
    """Convert covariance to correlation, clipping floating-point excursions."""
    variances = np.diag(covariance)
    if not np.isfinite(variances).all() or (variances <= 0.0).any():
        raise MicrobiomeSuiteError("SparCC covariance diagonal must be finite and positive.")
    scales = np.sqrt(variances)
    correlation = covariance / (scales[:, None] * scales[None, :])
    correlation = (correlation + correlation.T) / 2.0
    if not np.isfinite(correlation).all():
        raise MicrobiomeSuiteError("SparCC reconstructed correlation is not finite.")
    if (np.abs(correlation) > 1.0 + _CORRELATION_ROUNDOFF_TOLERANCE).any():
        raise MicrobiomeSuiteError("SparCC reconstructed correlation exceeds valid bounds.")
    correlation = np.clip(correlation, -1.0, 1.0)
    np.fill_diagonal(correlation, 1.0)
    return correlation


def _estimate_initial(compositions: np.ndarray) -> SparCCResult:
    """Compute one deterministic pre-exclusion SparCC estimate."""
    variation = _variation_matrix(_clr(compositions))
    coefficients, target = _basis_system(variation)
    variances = _solve_basis_variances(coefficients, target)
    covariance = _reconstruct_covariance(variation, variances)
    correlation = _covariance_to_correlation(covariance)
    return SparCCResult(covariance=covariance, correlation=correlation)
