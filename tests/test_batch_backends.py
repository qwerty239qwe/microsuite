from __future__ import annotations

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.backends import BATCH_BACKENDS, SUPPORTED_BACKENDS, resolve_backend


def test_every_backend_declares_a_known_scale() -> None:
    from microsuite.batch.value_type import VALUE_TYPES

    for backend in BATCH_BACKENDS.values():
        assert backend.value_type in VALUE_TYPES


def test_supported_backends_matches_the_table() -> None:
    assert set(SUPPORTED_BACKENDS) == set(BATCH_BACKENDS)
    assert "mmuphin" in SUPPORTED_BACKENDS


@pytest.mark.parametrize(
    ("backend", "value_type"),
    [
        ("mmuphin", "relative"),
        ("combat-seq", "counts"),
        ("conqur", "counts"),
        ("plsda-batch", "clr"),
        ("metadict", "relative"),
    ],
)
def test_declared_scales(backend: str, value_type: str) -> None:
    assert BATCH_BACKENDS[backend].value_type == value_type


def test_unknown_backend_lists_the_alternatives() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="mmuphin"):
        resolve_backend("combat", covariates=None, target=None)


def test_covariates_rejected_where_unsupported() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="covariates"):
        resolve_backend("plsda-batch", covariates=["sex"], target="disease")


def test_target_rejected_where_unsupported() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="target"):
        resolve_backend("mmuphin", covariates=None, target="disease")


def test_supervised_backend_requires_a_target() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="supervised"):
        resolve_backend("plsda-batch", covariates=None, target=None)


def test_supported_combination_resolves() -> None:
    backend = resolve_backend("mmuphin", covariates=["sex"], target=None)
    assert backend.name == "mmuphin"
    assert backend.script == "mmuphin"


def test_script_name_differs_from_backend_name_for_combat_seq() -> None:
    assert BATCH_BACKENDS["combat-seq"].script == "combat_seq"
    assert BATCH_BACKENDS["plsda-batch"].script == "plsda_batch"
