from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch import correct as correct_module
from microsuite.batch.correct import run_batch_correction


def _adata() -> ad.AnnData:
    adata = ad.AnnData(
        np.array(
            [
                [10.0, 30.0, 5.0],
                [20.0, 20.0, 5.0],
                [12.0, 28.0, 6.0],
                [18.0, 22.0, 4.0],
            ]
        )
    )
    adata.obs_names = ["s1", "s2", "s3", "s4"]
    adata.var_names = ["f1", "f2", "f3"]
    adata.obs = pd.DataFrame(
        {"run_id": ["A", "A", "B", "B"], "sex": ["m", "f", "m", "f"]},
        index=adata.obs_names,
    )
    return adata


def _confounded() -> ad.AnnData:
    adata = ad.AnnData(np.array([[10.0, 30.0, 5.0], [20.0, 20.0, 5.0]]))
    adata.obs_names = ["s1", "s2"]
    adata.var_names = ["f1", "f2", "f3"]
    adata.obs = pd.DataFrame({"run_id": ["A", "B"], "sex": ["m", "f"]}, index=adata.obs_names)
    return adata


def _fake_backend(write: dict[str, list[float]] | None = None, capture: dict | None = None):
    """Stand in for the R script: record the call, write a corrected table."""

    def _invoke(**kwargs) -> None:
        positional = kwargs["positional"]
        params = json.loads(Path(positional[2]).read_text(encoding="utf-8"))
        if capture is not None:
            capture.update(params=params, kwargs=kwargs)
        payload = write or {
            "f1": [11.0, 21.0, 13.0, 19.0],
            "f2": [31.0, 21.0, 29.0, 23.0],
            "f3": [6.0, 6.0, 7.0, 5.0],
        }
        n = len(next(iter(payload.values())))
        samples = [f"s{i + 1}" for i in range(n)]
        frame = pd.DataFrame(payload, index=samples).T
        frame.index.name = "feature_id"
        frame.to_csv(positional[3], sep="\t")

    return _invoke


def test_params_json_carries_batch_and_covariates(monkeypatch) -> None:
    capture: dict = {}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(capture=capture))
    run_batch_correction(_adata(), backend="mmuphin", batch="run_id", covariates=["sex"])
    assert capture["params"]["batch"] == "run_id"
    assert capture["params"]["covariates"] == ["sex"]
    assert capture["kwargs"]["script_name"] == "mmuphin"
    assert capture["kwargs"]["backend"] == "mmuphin"


def test_corrected_values_land_on_the_right_labels(monkeypatch) -> None:
    # The R script returns features as rows in an arbitrary order. A rebuild that
    # trusts position rather than labels puts f3's values under f1.
    payload = {
        "f3": [6.0, 6.0, 7.0, 5.0],
        "f1": [11.0, 21.0, 13.0, 19.0],
        "f2": [31.0, 21.0, 29.0, 23.0],
    }
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(write=payload))
    result = run_batch_correction(_adata(), backend="mmuphin", batch="run_id")
    assert list(result.var_names) == ["f1", "f2", "f3"]
    np.testing.assert_allclose(result.X[0], [11.0, 31.0, 6.0])
    np.testing.assert_allclose(result.X[1], [21.0, 21.0, 6.0])


def test_dropped_features_subset_var_rather_than_silently_realigning(monkeypatch) -> None:
    payload = {
        "f1": [11.0, 21.0, 13.0, 19.0],
        "f2": [31.0, 21.0, 29.0, 23.0],
    }
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(write=payload))
    result = run_batch_correction(_adata(), backend="mmuphin", batch="run_id")
    assert list(result.var_names) == ["f1", "f2"]
    assert result.shape == (4, 2)


def test_unknown_feature_in_output_raises(monkeypatch) -> None:
    payload = {
        "f1": [11.0, 21.0, 13.0, 19.0],
        "f9": [1.0, 1.0, 1.0, 1.0],
    }
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(write=payload))
    with pytest.raises(MicrobiomeSuiteError, match="f9"):
        run_batch_correction(_adata(), backend="mmuphin", batch="run_id")


def test_result_records_its_scale_and_provenance(monkeypatch) -> None:
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend())
    result = run_batch_correction(_adata(), backend="mmuphin", batch="run_id", covariates=["sex"])
    assert result.uns["microsuite"]["value_type"] == "relative"
    assert result.uns["microsuite"]["batch_correct"]["backend"] == "mmuphin"
    assert result.uns["microsuite"]["batch_correct"]["covariates"] == ["sex"]


def test_missing_batch_column_lists_available_columns() -> None:
    with pytest.raises(MicrobiomeSuiteError) as excinfo:
        run_batch_correction(_adata(), backend="mmuphin", batch="plate")
    assert "plate" in str(excinfo.value)
    assert "run_id" in str(excinfo.value)


def test_missing_covariate_column_raises() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="age"):
        run_batch_correction(_adata(), backend="mmuphin", batch="run_id", covariates=["age"])


def test_single_batch_level_raises() -> None:
    adata = _adata()
    adata.obs["run_id"] = ["A", "A", "A", "A"]
    with pytest.raises(MicrobiomeSuiteError, match="one batch"):
        run_batch_correction(adata, backend="mmuphin", batch="run_id")


def test_covariate_confounded_with_batch_raises() -> None:
    # 'sex' varies exactly with 'run_id' here, so no model can separate them.
    with pytest.raises(MicrobiomeSuiteError, match="confounded"):
        run_batch_correction(_confounded(), backend="mmuphin", batch="run_id", covariates=["sex"])


def test_combat_seq_declares_counts(monkeypatch) -> None:
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend())
    result = run_batch_correction(_adata(), backend="combat-seq", batch="run_id")
    assert result.uns["microsuite"]["value_type"] == "counts"


def test_combat_seq_uses_its_own_script_name(monkeypatch) -> None:
    capture: dict = {}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(capture=capture))
    run_batch_correction(_adata(), backend="combat-seq", batch="run_id")
    assert capture["kwargs"]["script_name"] == "combat_seq"
    assert capture["kwargs"]["backend"] == "combat-seq"


def test_plsda_batch_requires_a_target() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="supervised"):
        run_batch_correction(_adata(), backend="plsda-batch", batch="run_id")


def test_plsda_batch_rejects_covariates() -> None:
    with pytest.raises(MicrobiomeSuiteError, match="covariates"):
        run_batch_correction(
            _adata(),
            backend="plsda-batch",
            batch="run_id",
            covariates=["sex"],
            target="sex",
        )


def test_plsda_batch_declares_clr_and_passes_the_target(monkeypatch) -> None:
    capture: dict = {}
    monkeypatch.setattr(correct_module, "invoke_r_script", _fake_backend(capture=capture))
    result = run_batch_correction(_adata(), backend="plsda-batch", batch="run_id", target="sex")
    assert capture["params"]["target"] == "sex"
    assert result.uns["microsuite"]["value_type"] == "clr"
    assert result.uns["microsuite"]["batch_correct"]["target"] == "sex"
