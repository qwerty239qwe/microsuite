from __future__ import annotations

import importlib
from pathlib import Path

import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_tsv
from microsuite.methods.ml_longitudinal import (
    longitudinal,
    longitudinal_native,
    ml_classify,
    ml_classify_native,
)

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def longitudinal_adata():
    adata = fixture_adata().copy()
    adata.obs["subject"] = ["S1", "S1", "S2", "S2"]
    adata.obs["day"] = [0, 7, 0, 7]
    adata.obs["treatment"] = ["A", "A", "B", "B"]
    return adata


def fixture_table(tmp_path: Path) -> Path:
    table = tmp_path / "table.h5ad"
    write_h5ad(fixture_adata(), table)
    return table


def longitudinal_table(tmp_path: Path) -> Path:
    table = tmp_path / "longitudinal.h5ad"
    write_h5ad(longitudinal_adata(), table)
    return table


def test_randomforest_native_classifier_outputs_predictions_and_importance() -> None:
    predictions, importance = ml_classify_native(
        fixture_adata(),
        backend="randomforest",
        target="body_site",
        test_fraction=0.5,
        n_estimators=10,
        seed=2,
    )

    assert {"sample_id", "truth", "prediction", "correct", "confidence", "backend"}.issubset(
        predictions.columns
    )
    assert set(predictions["truth"]) == {"gut", "tongue"}
    assert {"feature", "importance", "backend"}.issubset(importance.columns)


def test_xgboost_missing_optional_dependency_reports_message(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_import_module = importlib.import_module

    def fake_import_module(name: str, package: str | None = None) -> object:
        if name == "xgboost":
            raise ImportError("missing")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", fake_import_module)

    with pytest.raises(MicrobiomeSuiteError, match="xgboost"):
        ml_classify_native(
            fixture_adata(),
            backend="xgboost",
            target="body_site",
            test_fraction=0.5,
        )


def test_ml_classify_writes_prediction_and_importance_files(tmp_path: Path) -> None:
    output = tmp_path / "predictions.tsv"
    importance = tmp_path / "importance.tsv"

    ml_classify(
        backend="randomforest",
        table=fixture_table(tmp_path),
        target="body_site",
        output=output,
        importance_output=importance,
        test_fraction=0.5,
        n_estimators=10,
        seed=1,
    )

    assert pd.read_csv(output, sep="\t").shape[0] == 2
    assert "importance" in pd.read_csv(importance, sep="\t").columns


def test_longitudinal_native_outputs_grouped_slopes() -> None:
    result = longitudinal_native(
        longitudinal_adata(),
        subject="subject",
        time="day",
        group="treatment",
        level="genus",
    )

    assert {"group", "feature", "n_subjects", "mean_slope"}.issubset(result.columns)
    assert sorted(result["group"].unique().tolist()) == ["A", "B"]
    assert (result["n_subjects"] == 1).all()


def test_longitudinal_cli_writes_slopes(tmp_path: Path) -> None:
    output = tmp_path / "slopes.tsv"
    result = CliRunner().invoke(
        app,
        [
            "ml",
            "longitudinal",
            "--backend",
            "native-time-series",
            "--table",
            str(longitudinal_table(tmp_path)),
            "--subject",
            "subject",
            "--time",
            "day",
            "--group",
            "treatment",
            "-o",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert pd.read_csv(output, sep="\t").shape[0] > 0


def test_ml_cli_classify_writes_predictions(tmp_path: Path) -> None:
    output = tmp_path / "predictions.tsv"
    result = CliRunner().invoke(
        app,
        [
            "ml",
            "classify",
            "--backend",
            "randomforest",
            "--table",
            str(fixture_table(tmp_path)),
            "--target",
            "body_site",
            "--test-fraction",
            "0.5",
            "-o",
            str(output),
        ],
    )

    assert result.exit_code == 0, result.stdout
    assert "prediction" in pd.read_csv(output, sep="\t").columns


def test_longitudinal_rejects_missing_metadata_column(tmp_path: Path) -> None:
    with pytest.raises(MicrobiomeSuiteError, match="missing"):
        longitudinal(
            backend="native-time-series",
            table=fixture_table(tmp_path),
            subject="missing",
            time="day",
            output=tmp_path / "slopes.tsv",
        )
