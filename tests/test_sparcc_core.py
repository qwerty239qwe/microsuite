from __future__ import annotations

from dataclasses import FrozenInstanceError
from inspect import signature

import numpy as np
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods._sparcc import (
    SparCCResult,
    _dirichlet_normalize,
    _validate_inputs,
    estimate_sparcc,
)


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

    with pytest.raises(NotImplementedError, match="basis algebra"):
        estimate_sparcc(valid_counts())
