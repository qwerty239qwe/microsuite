from __future__ import annotations

import anndata as ad
import numpy as np
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.value_type import (
    read_value_type,
    record_batch_correction,
    require_value_types,
)


def _adata() -> ad.AnnData:
    return ad.AnnData(np.array([[1.0, 2.0], [3.0, 4.0]]))


def test_absent_key_skips_the_check() -> None:
    # Any table written before 0.3.0 has no key at all. Behaviour must not change.
    require_value_types(_adata(), ("counts",), operation="rarefy")


def test_microsuite_key_without_value_type_skips_the_check() -> None:
    adata = _adata()
    adata.uns["microsuite"] = {"source": "tsv"}  # io/tsv.py writes this shape
    require_value_types(adata, ("counts",), operation="rarefy")


def test_matching_value_type_passes() -> None:
    adata = _adata()
    record_batch_correction(
        adata,
        value_type="counts",
        backend="combat-seq",
        batch="run_id",
        covariates=[],
        target=None,
    )
    require_value_types(adata, ("counts",), operation="rarefy")


def test_mismatched_value_type_names_the_producing_backend() -> None:
    adata = _adata()
    record_batch_correction(
        adata,
        value_type="clr",
        backend="plsda-batch",
        batch="run_id",
        covariates=[],
        target="disease",
    )
    with pytest.raises(MicrobiomeSuiteError) as excinfo:
        require_value_types(adata, ("counts",), operation="diff_abundance --backend ancombc")
    message = str(excinfo.value)
    assert "diff_abundance --backend ancombc" in message
    assert "clr" in message
    assert "plsda-batch" in message


def test_provenance_is_recorded_alongside_the_scale() -> None:
    adata = _adata()
    record_batch_correction(
        adata,
        value_type="relative",
        backend="mmuphin",
        batch="run_id",
        covariates=["sex", "age"],
        target=None,
    )
    assert read_value_type(adata) == "relative"
    provenance = adata.uns["microsuite"]["batch_correct"]
    assert provenance["backend"] == "mmuphin"
    assert provenance["batch"] == "run_id"
    assert provenance["covariates"] == ["sex", "age"]
    assert provenance["target"] is None


def test_recording_preserves_existing_microsuite_keys() -> None:
    adata = _adata()
    adata.uns["microsuite"] = {"source": "tsv"}
    record_batch_correction(
        adata,
        value_type="counts",
        backend="conqur",
        batch="run_id",
        covariates=[],
        target=None,
    )
    assert adata.uns["microsuite"]["source"] == "tsv"
    assert adata.uns["microsuite"]["value_type"] == "counts"


def test_unknown_value_type_is_rejected_at_write_time() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="abundance"):
        record_batch_correction(
            _adata(),
            value_type="abundance",
            backend="mmuphin",
            batch="run_id",
            covariates=[],
            target=None,
        )
