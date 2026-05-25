from __future__ import annotations

import importlib
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity._matrix import dense_counts
from microsuite.io.h5ad import read_h5ad

SUPPORTED_ML_BACKENDS = ("randomforest", "xgboost")
SUPPORTED_LONGITUDINAL_BACKENDS = ("native-time-series",)


def ml_classify(
    *,
    backend: str,
    table: Path,
    target: str,
    output: Path,
    importance_output: Path | None = None,
    test_fraction: float = 0.25,
    n_estimators: int = 100,
    seed: int = 0,
    force: bool = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    predictions, importance = ml_classify_native(
        adata,
        backend=backend,
        target=target,
        test_fraction=test_fraction,
        n_estimators=n_estimators,
        seed=seed,
    )
    predictions.to_csv(prepare_output(output, force=force), sep="\t", index=False)
    if importance_output is not None:
        importance.to_csv(prepare_output(importance_output, force=force), sep="\t", index=False)


def ml_classify_native(
    adata: ad.AnnData,
    *,
    backend: str,
    target: str,
    test_fraction: float = 0.25,
    n_estimators: int = 100,
    seed: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    backend = backend.lower()
    if backend not in SUPPORTED_ML_BACKENDS:
        raise MicrobiomeSuiteError(
            f"Unsupported ML backend '{backend}'. Choose one of: {', '.join(SUPPORTED_ML_BACKENDS)}"
        )
    if backend == "xgboost":
        return _classify_xgboost(
            adata,
            target=target,
            test_fraction=test_fraction,
            n_estimators=n_estimators,
            seed=seed,
        )
    return _classify_randomforest(
        adata,
        target=target,
        test_fraction=test_fraction,
        n_estimators=n_estimators,
        seed=seed,
    )


def longitudinal(
    *,
    backend: str,
    table: Path,
    subject: str,
    time: str,
    output: Path,
    group: str | None = None,
    level: str | None = None,
    force: bool = False,
) -> None:
    backend = backend.lower()
    if backend not in SUPPORTED_LONGITUDINAL_BACKENDS:
        raise MicrobiomeSuiteError(
            "Unsupported longitudinal backend "
            f"'{backend}'. Choose one of: {', '.join(SUPPORTED_LONGITUDINAL_BACKENDS)}"
        )
    result = longitudinal_native(
        read_h5ad(ensure_input(table)),
        subject=subject,
        time=time,
        group=group,
        level=level,
    )
    result.to_csv(prepare_output(output, force=force), sep="\t", index=False)


def longitudinal_native(
    adata: ad.AnnData,
    *,
    subject: str,
    time: str,
    group: str | None = None,
    level: str | None = None,
) -> pd.DataFrame:
    obs = pd.DataFrame(adata.obs)
    for column in [subject, time, *([group] if group is not None else [])]:
        if column not in obs.columns:
            raise MicrobiomeSuiteError(f"Sample metadata column not found: {column}")
    times = pd.to_numeric(obs[time], errors="raise").to_numpy(dtype=float)
    matrix, features = _feature_matrix(adata, level=level)
    rows = []
    if group is None:
        group_values = np.array(["all"] * adata.n_obs)
    else:
        group_values = obs[group].astype(str).to_numpy()
    for group_value in sorted(set(group_values)):
        group_mask = group_values == group_value
        for feature_index, feature in enumerate(features):
            subject_slopes = []
            for subject_value in sorted(set(obs.loc[group_mask, subject].astype(str))):
                mask = group_mask & (obs[subject].astype(str).to_numpy() == subject_value)
                if mask.sum() < 2:
                    continue
                slope = _slope(times[mask], matrix[mask, feature_index])
                subject_slopes.append(slope)
            if not subject_slopes:
                continue
            rows.append(
                {
                    "group": group_value,
                    "feature": feature,
                    "n_subjects": len(subject_slopes),
                    "mean_slope": float(np.mean(subject_slopes)),
                    "median_slope": float(np.median(subject_slopes)),
                    "min_slope": float(np.min(subject_slopes)),
                    "max_slope": float(np.max(subject_slopes)),
                }
            )
    return pd.DataFrame(rows).sort_values(
        ["group", "mean_slope", "feature"],
        ascending=[True, False, True],
    )


def _classify_randomforest(
    adata: ad.AnnData,
    *,
    target: str,
    test_fraction: float,
    n_estimators: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    x, y, sample_ids, feature_ids = _classification_inputs(adata, target=target)
    train_idx, test_idx = _split_indices(y, test_fraction=test_fraction, seed=seed)
    try:
        ensemble = importlib.import_module("sklearn.ensemble")
    except ImportError:
        predictions, importance = _nearest_centroid_classifier(
            x, y, sample_ids, feature_ids, train_idx=train_idx, test_idx=test_idx
        )
        predictions["backend"] = "randomforest-fallback"
        return predictions, importance

    model = ensemble.RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=seed,
        class_weight="balanced",
    )
    model.fit(x[train_idx], y[train_idx])
    predicted = model.predict(x[test_idx])
    probabilities = model.predict_proba(x[test_idx])
    confidence = probabilities.max(axis=1)
    predictions = _prediction_frame(sample_ids, y, test_idx, predicted, confidence, "randomforest")
    importance = pd.DataFrame(
        {
            "feature": feature_ids,
            "importance": model.feature_importances_,
            "backend": "randomforest",
        }
    ).sort_values(["importance", "feature"], ascending=[False, True])
    return predictions, importance.reset_index(drop=True)


def _classify_xgboost(
    adata: ad.AnnData,
    *,
    target: str,
    test_fraction: float,
    n_estimators: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    x, y, sample_ids, feature_ids = _classification_inputs(adata, target=target)
    train_idx, test_idx = _split_indices(y, test_fraction=test_fraction, seed=seed)
    try:
        xgboost = importlib.import_module("xgboost")
    except ImportError as exc:
        raise MicrobiomeSuiteError(
            "xgboost backend requires the optional 'xgboost' Python package."
        ) from exc
    labels = sorted(set(y))
    label_to_int = {label: i for i, label in enumerate(labels)}
    encoded_y = np.array([label_to_int[label] for label in y])
    model = xgboost.XGBClassifier(
        n_estimators=n_estimators,
        random_state=seed,
        eval_metric="mlogloss",
    )
    model.fit(x[train_idx], encoded_y[train_idx])
    predicted_encoded = model.predict(x[test_idx])
    predicted = np.array([labels[int(value)] for value in predicted_encoded])
    probabilities = model.predict_proba(x[test_idx])
    confidence = probabilities.max(axis=1)
    predictions = _prediction_frame(sample_ids, y, test_idx, predicted, confidence, "xgboost")
    importance = pd.DataFrame(
        {
            "feature": feature_ids,
            "importance": model.feature_importances_,
            "backend": "xgboost",
        }
    ).sort_values(["importance", "feature"], ascending=[False, True])
    return predictions, importance.reset_index(drop=True)


def _classification_inputs(
    adata: ad.AnnData,
    *,
    target: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    obs = pd.DataFrame(adata.obs)
    if target not in obs.columns:
        raise MicrobiomeSuiteError(f"Sample metadata target not found: {target}")
    y = obs[target].astype(str).to_numpy()
    if len(set(y)) < 2:
        raise MicrobiomeSuiteError("Classification target must contain at least two classes.")
    counts = dense_counts(adata).astype(float)
    totals = counts.sum(axis=1, keepdims=True)
    x = np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0)
    return x, y, adata.obs_names.astype(str).to_numpy(), adata.var_names.astype(str).tolist()


def _split_indices(
    y: np.ndarray, *, test_fraction: float, seed: int
) -> tuple[np.ndarray, np.ndarray]:
    if not 0 < test_fraction < 1:
        raise MicrobiomeSuiteError("--test-fraction must be greater than 0 and less than 1.")
    rng = np.random.default_rng(seed)
    train = []
    test = []
    for label in sorted(set(y)):
        indices = np.flatnonzero(y == label)
        rng.shuffle(indices)
        if len(indices) < 2:
            raise MicrobiomeSuiteError("Each class must have at least two samples for splitting.")
        n_test = max(1, int(round(len(indices) * test_fraction)))
        n_test = min(n_test, len(indices) - 1)
        test.extend(indices[:n_test].tolist())
        train.extend(indices[n_test:].tolist())
    return np.array(sorted(train)), np.array(sorted(test))


def _nearest_centroid_classifier(
    x: np.ndarray,
    y: np.ndarray,
    sample_ids: np.ndarray,
    feature_ids: list[str],
    *,
    train_idx: np.ndarray,
    test_idx: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    labels = sorted(set(y[train_idx]))
    centroids = np.vstack([x[train_idx][y[train_idx] == label].mean(axis=0) for label in labels])
    distances = np.vstack([np.linalg.norm(centroids - row, axis=1) for row in x[test_idx]])
    nearest = distances.argmin(axis=1)
    predicted = np.array([labels[index] for index in nearest])
    confidence = 1.0 / (1.0 + distances.min(axis=1))
    predictions = _prediction_frame(
        sample_ids,
        y,
        test_idx,
        predicted,
        confidence,
        "randomforest-fallback",
    )
    spread = np.var(centroids, axis=0)
    importance = pd.DataFrame(
        {
            "feature": feature_ids,
            "importance": spread,
            "backend": "randomforest-fallback",
        }
    ).sort_values(["importance", "feature"], ascending=[False, True])
    return predictions, importance.reset_index(drop=True)


def _prediction_frame(
    sample_ids: np.ndarray,
    y: np.ndarray,
    test_idx: np.ndarray,
    predicted: np.ndarray,
    confidence: np.ndarray,
    backend: str,
) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "sample_id": sample_ids[test_idx],
            "truth": y[test_idx],
            "prediction": predicted,
            "correct": predicted == y[test_idx],
            "confidence": confidence,
            "backend": backend,
        }
    )


def _feature_matrix(adata: ad.AnnData, *, level: str | None) -> tuple[np.ndarray, list[str]]:
    counts = dense_counts(adata).astype(float)
    totals = counts.sum(axis=1, keepdims=True)
    relative = np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0)
    if level is None:
        return relative, adata.var_names.astype(str).tolist()
    var = pd.DataFrame(adata.var)
    if level not in var.columns:
        raise MicrobiomeSuiteError(f"Taxonomy level not found in table: {level}")
    taxa = var[level].astype(str).replace("", "Unclassified")
    frame = pd.DataFrame(relative, index=adata.obs_names.astype(str), columns=taxa)
    grouped = frame.T.groupby(level=0).sum().T
    return grouped.to_numpy(dtype=float), grouped.columns.astype(str).tolist()


def _slope(x: np.ndarray, y: np.ndarray) -> float:
    x_centered = x - x.mean()
    denominator = np.square(x_centered).sum()
    if denominator == 0:
        return 0.0
    return float((x_centered * (y - y.mean())).sum() / denominator)
