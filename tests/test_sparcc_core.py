from __future__ import annotations

from dataclasses import FrozenInstanceError
from inspect import signature
from pathlib import Path

import numpy as np
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods._sparcc import (
    SparCCResult,
    _basis_system,
    _clr,
    _covariance_to_correlation,
    _dirichlet_normalize,
    _estimate_initial,
    _reconstruct_covariance,
    _solve_basis_variances,
    _validate_inputs,
    _variation_matrix,
    estimate_sparcc,
)

SPARCC_FIXTURES = Path(__file__).parent / "fixtures" / "sparcc"


def valid_counts() -> np.ndarray:
    return np.array([[4, 0, 2], [1, 3, 1], [0, 2, 5]], dtype=np.int64)


def validate(counts: np.ndarray, **overrides: object) -> np.ndarray:
    parameters: dict[str, object] = {
        "iterations": 20,
        "inner_iterations": 10,
        "exclusion_threshold": 0.1,
        "pseudocount": 1.0,
        "seed": 0,
    }
    parameters.update(overrides)
    return _validate_inputs(counts, **parameters)  # type: ignore[arg-type]


def test_sparcc_result_is_frozen_and_estimator_signature_is_stable() -> None:
    covariance = np.eye(3)
    result = SparCCResult(covariance=covariance, correlation=np.eye(3))

    with pytest.raises(FrozenInstanceError):
        result.covariance = np.zeros((3, 3))

    parameters = signature(estimate_sparcc).parameters
    assert list(parameters) == [
        "counts",
        "iterations",
        "inner_iterations",
        "exclusion_threshold",
        "pseudocount",
        "seed",
    ]
    assert [parameter.default for parameter in list(parameters.values())[1:]] == [
        20,
        10,
        0.1,
        1.0,
        0,
    ]


@pytest.mark.parametrize(
    "counts",
    [
        np.array([1, 2, 3]),
        np.ones((1, 3)),
        np.ones((2, 2)),
        np.array([[1, 2, np.nan], [2, 1, 3]]),
        np.array([[1, 2, np.inf], [2, 1, 3]]),
        np.array([[1, 2, -1], [2, 1, 3]]),
        np.array([[1, 2, 1.5], [2, 1, 3]]),
        np.array([[0, 0, 0], [2, 1, 3]]),
        np.array([[1, 0, 2], [2, 0, 3]]),
        np.array([["one", "2", "3"], ["2", "1", "3"]]),
        np.array([[True, False, True], [False, True, True]]),
        np.array([[1 + 1j, 2, 3], [2, 1, 3]]),
    ],
)
def test_invalid_count_matrices_raise_project_error(counts: np.ndarray) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        validate(counts)


@pytest.mark.parametrize("name", ["iterations", "inner_iterations"])
@pytest.mark.parametrize("value", [0, -1, 1.0, True])
def test_iteration_parameters_require_positive_integers(name: str, value: object) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        validate(valid_counts(), **{name: value})


@pytest.mark.parametrize("value", [-0.01, 1.01, np.nan, np.inf, True, "0.1"])
def test_exclusion_threshold_must_be_finite_and_in_range(value: object) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        validate(valid_counts(), exclusion_threshold=value)


@pytest.mark.parametrize("value", [0.0, -1.0, np.nan, np.inf, True, "1.0"])
def test_pseudocount_must_be_positive_and_finite(value: object) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        validate(valid_counts(), pseudocount=value)


@pytest.mark.parametrize("value", [-1, 1.0, True, "1"])
def test_seed_must_be_a_nonnegative_integer(value: object) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        validate(valid_counts(), seed=value)


def test_exact_integer_floats_are_accepted_without_mutating_caller() -> None:
    counts = valid_counts().astype(np.float64)
    before = counts.copy()

    validated = validate(counts)

    np.testing.assert_array_equal(counts, before)
    np.testing.assert_array_equal(validated, counts)
    assert validated.dtype == np.float64
    assert not np.shares_memory(validated, counts)

    with pytest.raises(NotImplementedError):
        estimate_sparcc(counts)
    np.testing.assert_array_equal(counts, before)


def draw(counts: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return _dirichlet_normalize(counts.astype(np.float64), pseudocount=1.0, rng=rng)


def test_dirichlet_draws_are_positive_and_row_normalized() -> None:
    compositions = draw(valid_counts(), seed=12)

    assert (compositions > 0.0).all()
    np.testing.assert_allclose(compositions.sum(axis=1), 1.0, rtol=0.0, atol=1e-15)


def test_dirichlet_draws_repeat_for_same_seed_and_change_for_different_seed() -> None:
    first = draw(valid_counts(), seed=12)
    repeated = draw(valid_counts(), seed=12)
    different = draw(valid_counts(), seed=13)

    np.testing.assert_array_equal(first, repeated)
    assert not np.array_equal(first, different)


def test_dirichlet_draws_do_not_touch_numpy_global_rng() -> None:
    np.random.seed(8675309)
    expected = np.random.random(5)
    np.random.seed(8675309)

    draw(valid_counts(), seed=12)

    np.testing.assert_array_equal(np.random.random(5), expected)


def test_estimator_validates_before_the_unimplemented_algebra() -> None:
    with pytest.raises(MicrobiomeSuiteError):
        estimate_sparcc(np.ones((1, 3)))

    with pytest.raises(NotImplementedError, match="iterative exclusion"):
        estimate_sparcc(valid_counts())


def test_clr_is_centered_within_each_sample() -> None:
    compositions = np.array([[0.1, 0.2, 0.7], [0.6, 0.3, 0.1]])

    clr_values = _clr(compositions)

    np.testing.assert_allclose(clr_values.mean(axis=1), 0.0, rtol=0.0, atol=1e-15)


def test_variation_uses_sample_variance_and_has_aitchison_invariants() -> None:
    compositions = np.array([[0.1, 0.2, 0.7], [0.6, 0.3, 0.1], [0.25, 0.25, 0.5]])
    logged = np.log(compositions)
    expected = np.empty((3, 3))
    population = np.empty((3, 3))
    for left in range(3):
        for right in range(3):
            log_ratio = logged[:, left] - logged[:, right]
            expected[left, right] = np.var(log_ratio, ddof=1)
            population[left, right] = np.var(log_ratio, ddof=0)

    variation = _variation_matrix(_clr(compositions))

    np.testing.assert_allclose(variation, expected, rtol=1e-14, atol=1e-15)
    assert not np.allclose(variation, population)
    np.testing.assert_allclose(variation, variation.T, rtol=0.0, atol=0.0)
    np.testing.assert_array_equal(np.diag(variation), np.zeros(3))
    assert variation.min() >= -1e-15


def test_initial_basis_system_is_hand_computable_and_floors_variances() -> None:
    variation = np.array([[0.0, 3.0, 4.0], [3.0, 0.0, 5.0], [4.0, 5.0, 0.0]])

    coefficients, target = _basis_system(variation)
    variances = _solve_basis_variances(coefficients, target)

    np.testing.assert_array_equal(
        coefficients, np.array([[2.0, 1.0, 1.0], [1.0, 2.0, 1.0], [1.0, 1.0, 2.0]])
    )
    np.testing.assert_array_equal(target, np.array([7.0, 8.0, 9.0]))
    np.testing.assert_allclose(variances, np.array([1.0, 2.0, 3.0]))
    np.testing.assert_array_equal(
        _solve_basis_variances(np.eye(3), np.array([-1.0, 0.0, 2.0])),
        np.array([1e-4, 1e-4, 2.0]),
    )


def test_covariance_reconstruction_is_symmetric_bounded_and_has_expected_diagonal() -> None:
    variances = np.array([1.0, 4.0, 9.0])
    variation = np.array([[0.0, 3.0, 12.0], [3.0, 0.0, 7.0], [12.0, 7.0, 0.0]])
    expected_covariance = np.array([[1.0, 1.0, -1.0], [1.0, 4.0, 3.0], [-1.0, 3.0, 9.0]])

    covariance = _reconstruct_covariance(variation, variances)
    correlation = _covariance_to_correlation(covariance)

    np.testing.assert_array_equal(covariance, expected_covariance)
    np.testing.assert_array_equal(covariance, covariance.T)
    np.testing.assert_array_equal(np.diag(covariance), variances)
    np.testing.assert_array_equal(np.diag(correlation), np.ones(3))
    assert np.max(np.abs(correlation)) <= 1.0
    roundoff_covariance = np.array([[1.0, 1.0 + 5e-13], [1.0 + 5e-13, 1.0]])
    np.testing.assert_array_equal(_covariance_to_correlation(roundoff_covariance), np.ones((2, 2)))


def test_singular_basis_system_uses_generalized_inverse(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    original_pinv = np.linalg.pinv
    solve_calls = 0
    pinv_calls = 0

    def recording_solve(coefficients: np.ndarray, target: np.ndarray) -> np.ndarray:
        nonlocal solve_calls
        solve_calls += 1
        raise np.linalg.LinAlgError("singular test system")

    def recording_pinv(coefficients: np.ndarray) -> np.ndarray:
        nonlocal pinv_calls
        pinv_calls += 1
        return original_pinv(coefficients)

    monkeypatch.setattr(np.linalg, "solve", recording_solve)
    monkeypatch.setattr(np.linalg, "pinv", recording_pinv)
    coefficients = np.ones((3, 3))
    variances = _solve_basis_variances(coefficients, np.ones(3))
    covariance = _reconstruct_covariance(np.zeros((3, 3)), variances)
    correlation = _covariance_to_correlation(covariance)

    assert solve_calls == 1
    assert pinv_calls == 1
    assert np.isfinite(correlation).all()
    np.testing.assert_array_equal(correlation, correlation.T)


def test_initial_estimate_matches_pinned_spieceasi_pre_exclusion_result() -> None:
    compositions = np.loadtxt(
        SPARCC_FIXTURES / "inner_compositions.tsv",
        delimiter="\t",
        skiprows=1,
        usecols=range(1, 7),
    )
    expected = np.loadtxt(
        SPARCC_FIXTURES / "inner_initial_reference_cor.tsv",
        delimiter="\t",
        skiprows=1,
        usecols=range(1, 7),
    )

    result = _estimate_initial(compositions)

    np.testing.assert_allclose(result.correlation, expected, rtol=1e-10, atol=1e-12)
