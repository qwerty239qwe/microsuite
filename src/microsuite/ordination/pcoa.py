from __future__ import annotations

import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError


def pcoa(distance_matrix: pd.DataFrame, *, dimensions: int = 3) -> pd.DataFrame:
    if distance_matrix.shape[0] != distance_matrix.shape[1]:
        raise MicrobiomeSuiteError("Distance matrix must be square.")
    if list(distance_matrix.index.astype(str)) != list(distance_matrix.columns.astype(str)):
        raise MicrobiomeSuiteError("Distance matrix row and column IDs must match in order.")

    distances = distance_matrix.to_numpy(dtype=np.float64)
    if not np.allclose(distances, distances.T):
        raise MicrobiomeSuiteError("Distance matrix must be symmetric.")

    n_samples = distances.shape[0]
    centering = np.eye(n_samples) - np.ones((n_samples, n_samples)) / n_samples
    gram = -0.5 * centering @ np.square(distances) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = np.maximum(eigenvalues[:dimensions], 0.0)
    coordinates = eigenvectors[:, :dimensions] * np.sqrt(positive)

    rows: dict[str, object] = {"sample_id": distance_matrix.index.astype(str)}
    for index in range(dimensions):
        rows[f"PC{index + 1}"] = coordinates[:, index]

    total_positive = np.maximum(eigenvalues, 0.0).sum()
    for index in range(dimensions):
        rows[f"PC{index + 1}_variance"] = (
            0.0 if total_positive == 0 else positive[index] / total_positive
        )
    return pd.DataFrame(rows)
