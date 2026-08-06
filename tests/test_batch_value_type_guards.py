from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.value_type import record_batch_correction
from microsuite.io.h5ad import write_h5ad
from microsuite.methods.normalize import normalize_native
from microsuite.methods.rarefy import rarefy_native


def _table(value_type: str | None, backend: str = "mmuphin") -> ad.AnnData:
    adata = ad.AnnData(np.array([[10.0, 30.0], [20.0, 20.0]]))
    adata.obs_names = ["s1", "s2"]
    adata.var_names = ["f1", "f2"]
    if value_type is not None:
        record_batch_correction(
            adata,
            value_type=value_type,
            backend=backend,
            batch="run_id",
            covariates=[],
            target=None,
        )
    return adata


@pytest.mark.parametrize("value_type", ["relative", "clr"])
def test_rarefy_rejects_non_counts(value_type: str) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="rarefy"):
        rarefy_native(_table(value_type), depth=10)


def test_rarefy_accepts_counts_and_unmarked_tables() -> None:
    rarefy_native(_table("counts"), depth=10)
    rarefy_native(_table(None), depth=10)


@pytest.mark.parametrize("method", ["relative", "total-sum"])
@pytest.mark.parametrize("value_type", ["relative", "clr"])
def test_normalize_rejects_already_scaled_input(method: str, value_type: str) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="normalize"):
        normalize_native(_table(value_type), method=method)


def test_normalize_clr_rejects_clr_but_accepts_relative() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="normalize"):
        normalize_native(_table("clr"), method="clr")
    normalize_native(_table("relative"), method="clr")


def test_prevalence_filter_is_scale_agnostic() -> None:
    # Filtering by prevalence is valid at any scale, so it carries no guard.
    normalize_native(_table("clr"), method="prevalence-filter")


@pytest.mark.parametrize("backend", ["ancombc", "aldex2"])
def test_diff_abundance_count_backends_reject_clr(
    backend: str, tmp_path: Path, monkeypatch
) -> None:
    from microsuite.methods import diff_abundance as module

    table = tmp_path / "corrected.h5ad"
    write_h5ad(_table("clr", backend="plsda-batch"), table)
    monkeypatch.setattr(
        module, "run_ancombc", lambda *a, **kw: pytest.fail("backend must not be invoked")
    )
    monkeypatch.setattr(
        module,
        "run_r_diffab_backend",
        lambda *a, **kw: pytest.fail("backend must not be invoked"),
    )
    with pytest.raises(MicrobiomeSuiteError, match="plsda-batch"):
        module.diff_abundance(backend=backend, table=table, group="g", output=tmp_path / "out.tsv")


@pytest.mark.parametrize("backend", ["maaslin2", "lefse"])
def test_diff_abundance_internally_normalizing_backends_do_not_check(
    backend: str, tmp_path: Path, monkeypatch
) -> None:
    from microsuite.methods import diff_abundance as module

    invoked: dict = {}
    table = tmp_path / "corrected.h5ad"
    write_h5ad(_table("clr", backend="plsda-batch"), table)
    monkeypatch.setattr(module, "run_r_diffab_backend", lambda *a, **kw: invoked.update(ran=True))
    module.diff_abundance(backend=backend, table=table, group="g", output=tmp_path / "out.tsv")
    assert invoked["ran"] is True
