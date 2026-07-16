from __future__ import annotations

from dataclasses import dataclass

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts
from microsuite.diversity.alpha import _compute as _compute_alpha
from microsuite.ordination.pcoa import pcoa


@dataclass(frozen=True)
class PermutationResult:
    statistic: float
    p_value: float
    permutations: int


def beta_significance(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    column: str,
    method: str = "permanova",
    permutations: int = 999,
    seed: int = 0,
) -> pd.DataFrame:
    samples, distances, labels = _aligned_distance_and_labels(distance_matrix, metadata, column)
    method = method.lower()
    if method == "permanova":
        result = _permanova(distances, labels, permutations=permutations, seed=seed)
        statistic_name = "pseudo_f"
    elif method == "permdisp":
        result = _permdisp(distances, labels, permutations=permutations, seed=seed)
        statistic_name = "f_value"
    elif method == "anosim":
        result = _anosim(distances, labels, permutations=permutations, seed=seed)
        statistic_name = "r"
    else:
        raise MicrobiomeSuiteError("--method must be permanova, permdisp, or anosim.")
    return pd.DataFrame(
        [
            {
                "method": method,
                "metadata_column": column,
                "n_samples": len(samples),
                "n_groups": len(set(labels)),
                statistic_name: result.statistic,
                "p_value": result.p_value,
                "permutations": result.permutations,
                "permutation_scheme": "unrestricted",
            }
        ]
    )


def mantel_test(
    matrix_a: pd.DataFrame,
    matrix_b: pd.DataFrame,
    *,
    method: str = "pearson",
    permutations: int = 999,
    seed: int = 0,
) -> pd.DataFrame:
    matrix_b_samples = set(matrix_b.index.astype(str))
    samples = [sample for sample in matrix_a.index.astype(str) if sample in matrix_b_samples]
    if len(samples) < 3:
        raise MicrobiomeSuiteError("Mantel test requires at least three shared samples.")
    a = _distance_array(matrix_a.loc[samples, samples])
    b = _distance_array(matrix_b.loc[samples, samples])
    observed = _vector_correlation(_lower_triangle(a), _lower_triangle(b), method=method)
    rng = np.random.default_rng(seed)
    exceedances = 1
    for _ in range(permutations):
        permuted = rng.permutation(a.shape[0])
        statistic = _vector_correlation(
            _lower_triangle(a[permuted][:, permuted]),
            _lower_triangle(b),
            method=method,
        )
        if abs(statistic) >= abs(observed):
            exceedances += 1
    return pd.DataFrame(
        [
            {
                "method": f"mantel-{method.lower()}",
                "n_samples": len(samples),
                "r": observed,
                "p_value": exceedances / (permutations + 1),
                "permutations": permutations,
            }
        ]
    )


def gamma_diversity(
    adata: ad.AnnData,
    *,
    group: str,
    metric: str = "observed_features",
) -> pd.DataFrame:
    metric = metric.lower().replace("-", "_")
    obs = pd.DataFrame(adata.obs)
    if group not in obs.columns:
        raise MicrobiomeSuiteError(f"Sample metadata group not found: {group}")
    counts = dense_counts(adata)
    rows = []
    for value in sorted(obs[group].astype(str).unique()):
        mask = obs[group].astype(str).to_numpy() == value
        pooled = counts[mask].sum(axis=0, keepdims=True)
        rows.append(
            {
                "group": value,
                "n_samples": int(mask.sum()),
                metric: float(_compute_alpha(pooled, metric)[0]),
            }
        )
    return pd.DataFrame(rows)


def beta_turnover(
    adata: ad.AnnData,
    *,
    level: str | None = None,
) -> pd.DataFrame:
    labels, present = _presence_by_label(adata, level=level)
    rows = []
    for i, label_a in enumerate(labels):
        for j in range(i + 1, len(labels)):
            label_b = labels[j]
            shared = int(np.logical_and(present[i], present[j]).sum())
            only_a = int(np.logical_and(present[i], ~present[j]).sum())
            only_b = int(np.logical_and(~present[i], present[j]).sum())
            denom = (2 * shared) + only_a + only_b
            sorensen = 0.0 if denom == 0 else (only_a + only_b) / denom
            simpson = (
                0.0
                if shared + min(only_a, only_b) == 0
                else min(only_a, only_b) / (shared + min(only_a, only_b))
            )
            rows.append(
                {
                    "sample_a": label_a,
                    "sample_b": label_b,
                    "shared_taxa": shared,
                    "unique_to_a": only_a,
                    "unique_to_b": only_b,
                    "sorensen_dissimilarity": sorensen,
                    "turnover_component": simpson,
                    "nestedness_component": sorensen - simpson,
                }
            )
    return pd.DataFrame(rows)


def taxa_turnover(
    adata: ad.AnnData,
    *,
    group: str,
    level: str | None = None,
) -> pd.DataFrame:
    obs = pd.DataFrame(adata.obs)
    if group not in obs.columns:
        raise MicrobiomeSuiteError(f"Sample metadata group not found: {group}")
    labels, present = _presence_by_label(adata, level=level)
    groups = obs[group].astype(str).to_numpy()
    group_labels = sorted(set(groups))
    group_presence = np.vstack([present[groups == value].any(axis=0) for value in group_labels])
    rows = []
    for i, group_a in enumerate(group_labels):
        for j in range(i + 1, len(group_labels)):
            group_b = group_labels[j]
            shared = int(np.logical_and(group_presence[i], group_presence[j]).sum())
            only_a = int(np.logical_and(group_presence[i], ~group_presence[j]).sum())
            only_b = int(np.logical_and(~group_presence[i], group_presence[j]).sum())
            union = shared + only_a + only_b
            rows.append(
                {
                    "group_a": group_a,
                    "group_b": group_b,
                    "shared_taxa": shared,
                    "unique_to_a": only_a,
                    "unique_to_b": only_b,
                    "jaccard_turnover": 0.0 if union == 0 else (only_a + only_b) / union,
                    "taxa": ",".join(labels),
                }
            )
    return pd.DataFrame(rows)


def constrained_ordination(
    adata: ad.AnnData,
    *,
    constraints: list[str],
    method: str = "rda",
    dimensions: int = 2,
) -> pd.DataFrame:
    if not constraints:
        raise MicrobiomeSuiteError("At least one constraint column is required.")
    obs = pd.DataFrame(adata.obs)
    missing = [column for column in constraints if column not in obs.columns]
    if missing:
        raise MicrobiomeSuiteError(f"Constraint column not found: {', '.join(missing)}")
    method = method.lower()
    if method == "rda":
        response = _hellinger(dense_counts(adata))
    elif method == "cca":
        response = _chi_square_standardize(dense_counts(adata))
    elif method == "db-rda":
        response = pcoa(_bray_curtis_frame(adata), dimensions=min(adata.n_obs - 1, adata.n_vars))
        response = response.filter(regex=r"^PC").to_numpy(dtype=float)
    else:
        raise MicrobiomeSuiteError("--method must be rda, db-rda, or cca.")
    design = _design_matrix(obs, constraints)
    fitted = _fit_multivariate_response(design, response)
    coords, variance = _principal_coordinates(fitted, dimensions=dimensions)
    result = pd.DataFrame({"sample_id": adata.obs_names.astype(str)})
    for i in range(coords.shape[1]):
        result[f"Axis{i + 1}"] = coords[:, i]
        result[f"Axis{i + 1}_variance"] = variance[i]
    result["method"] = method
    return result


def _aligned_distance_and_labels(
    distance_matrix: pd.DataFrame, metadata: pd.DataFrame, column: str
) -> tuple[list[str], np.ndarray, np.ndarray]:
    if column not in metadata.columns:
        raise MicrobiomeSuiteError(f"Metadata column not found: {column}")
    matrix = distance_matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    meta = metadata.copy()
    meta.index = meta.index.astype(str)
    samples = [sample for sample in matrix.index if sample in set(meta.index)]
    if len(samples) < 3:
        raise MicrobiomeSuiteError("At least three samples with metadata are required.")
    labels = meta.loc[samples, column].astype(str).to_numpy()
    if len(set(labels)) < 2:
        raise MicrobiomeSuiteError("At least two groups are required.")
    return samples, _distance_array(matrix.loc[samples, samples]), labels


def _permanova(
    distances: np.ndarray, labels: np.ndarray, *, permutations: int, seed: int
) -> PermutationResult:
    observed = _permanova_statistic(distances, labels)
    rng = np.random.default_rng(seed)
    exceedances = 1
    for _ in range(permutations):
        if _permanova_statistic(distances, rng.permutation(labels)) >= observed:
            exceedances += 1
    return PermutationResult(observed, exceedances / (permutations + 1), permutations)


def _permanova_statistic(distances: np.ndarray, labels: np.ndarray) -> float:
    n = distances.shape[0]
    groups = sorted(set(labels))
    squared = np.square(distances)
    total_ss = squared.sum() / (2 * n)
    within_ss = 0.0
    for group in groups:
        indices = np.flatnonzero(labels == group)
        within_ss += squared[np.ix_(indices, indices)].sum() / (2 * len(indices))
    between_ss = total_ss - within_ss
    df_between = len(groups) - 1
    df_within = n - len(groups)
    if df_within <= 0 or within_ss <= 0:
        return 0.0
    return (between_ss / df_between) / (within_ss / df_within)


def _permdisp(
    distances: np.ndarray, labels: np.ndarray, *, permutations: int, seed: int
) -> PermutationResult:
    coordinates = _positive_pcoa_coordinates(distances)
    observed = _permdisp_statistic(coordinates, labels)
    rng = np.random.default_rng(seed)
    exceedances = 1
    for _ in range(permutations):
        if _permdisp_statistic(coordinates, rng.permutation(labels)) >= observed:
            exceedances += 1
    return PermutationResult(observed, exceedances / (permutations + 1), permutations)


def _positive_pcoa_coordinates(distances: np.ndarray) -> np.ndarray:
    n = distances.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gram = -0.5 * centering @ np.square(distances) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    positive = eigenvalues > max(float(np.max(np.abs(eigenvalues))) * 1e-12, 0.0)
    if not np.any(positive):
        return np.zeros((n, 1), dtype=float)
    return eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])


def _permdisp_statistic(coordinates: np.ndarray, labels: np.ndarray) -> float:
    groups = sorted(set(labels))
    dispersions = np.zeros(coordinates.shape[0], dtype=float)
    for group in groups:
        indices = np.flatnonzero(labels == group)
        centroid = coordinates[indices].mean(axis=0)
        dispersions[indices] = np.linalg.norm(coordinates[indices] - centroid, axis=1)
    grand_mean = float(dispersions.mean())
    between_ss = 0.0
    within_ss = 0.0
    for group in groups:
        values = dispersions[labels == group]
        group_mean = float(values.mean())
        between_ss += len(values) * (group_mean - grand_mean) ** 2
        within_ss += float(np.square(values - group_mean).sum())
    df_between = len(groups) - 1
    df_within = len(dispersions) - len(groups)
    if df_between <= 0 or df_within <= 0:
        return 0.0
    if within_ss <= np.finfo(float).eps:
        return float("inf") if between_ss > np.finfo(float).eps else 0.0
    return (between_ss / df_between) / (within_ss / df_within)


def _anosim(
    distances: np.ndarray, labels: np.ndarray, *, permutations: int, seed: int
) -> PermutationResult:
    observed = _anosim_statistic(distances, labels)
    rng = np.random.default_rng(seed)
    exceedances = 1
    for _ in range(permutations):
        if _anosim_statistic(distances, rng.permutation(labels)) >= observed:
            exceedances += 1
    return PermutationResult(observed, exceedances / (permutations + 1), permutations)


def _anosim_statistic(distances: np.ndarray, labels: np.ndarray) -> float:
    n = distances.shape[0]
    lower = np.tril_indices(n, k=-1)
    ranks = pd.Series(distances[lower]).rank(method="average").to_numpy()
    same = labels[lower[0]] == labels[lower[1]]
    if same.all() or (~same).all():
        return 0.0
    return float((ranks[~same].mean() - ranks[same].mean()) / (n * (n - 1) / 4))


def _presence_by_label(adata: ad.AnnData, *, level: str | None) -> tuple[list[str], np.ndarray]:
    counts = dense_counts(adata)
    labels = adata.obs_names.astype(str).tolist()
    if level is None:
        return labels, counts > 0
    var = pd.DataFrame(adata.var)
    if level not in var.columns:
        raise MicrobiomeSuiteError(f"Taxonomy level not found in table: {level}")
    taxa = var[level].astype(str).replace("", "Unclassified")
    frame = pd.DataFrame(counts, index=labels, columns=taxa)
    grouped = frame.T.groupby(level=0).sum().T
    return labels, grouped.to_numpy() > 0


def _design_matrix(obs: pd.DataFrame, constraints: list[str]) -> np.ndarray:
    parts = []
    for column in constraints:
        values = obs[column]
        if pd.api.types.is_numeric_dtype(values):
            centered = values.astype(float).to_numpy()
            centered = centered - centered.mean()
            parts.append(pd.DataFrame({column: centered}, index=obs.index))
        else:
            encoded = pd.get_dummies(
                values.astype(str),
                prefix=column,
                drop_first=True,
                dtype=float,
            )
            parts.append(encoded)
    if not parts:
        raise MicrobiomeSuiteError("No usable constraint columns were provided.")
    matrix = pd.concat(parts, axis=1).to_numpy(dtype=float)
    return np.column_stack([np.ones(matrix.shape[0]), matrix])


def _fit_multivariate_response(design: np.ndarray, response: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(design, response, rcond=None)
    return design @ beta


def _principal_coordinates(matrix: np.ndarray, *, dimensions: int) -> tuple[np.ndarray, np.ndarray]:
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    _, singular_values, vt = np.linalg.svd(centered, full_matrices=False)
    dims = max(1, min(dimensions, vt.shape[0]))
    coords = centered @ vt[:dims].T
    eigenvalues = np.square(singular_values[:dims])
    total = np.square(singular_values).sum()
    variance = np.zeros(dims) if total == 0 else eigenvalues / total
    return coords, variance


def _hellinger(counts: np.ndarray) -> np.ndarray:
    totals = counts.sum(axis=1, keepdims=True)
    proportions = np.divide(
        counts,
        totals,
        out=np.zeros_like(counts, dtype=float),
        where=totals > 0,
    )
    return np.sqrt(proportions)


def _chi_square_standardize(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    if total == 0:
        return counts.astype(float)
    proportions = counts / total
    row_weights = proportions.sum(axis=1, keepdims=True)
    column_weights = proportions.sum(axis=0, keepdims=True)
    expected = row_weights @ column_weights
    standardized = np.divide(
        proportions - expected,
        np.sqrt(expected),
        out=np.zeros_like(proportions, dtype=float),
        where=expected > 0,
    )
    return standardized


def _bray_curtis_frame(adata: ad.AnnData) -> pd.DataFrame:
    from microsuite.diversity.beta import beta_diversity

    return beta_diversity(adata, "bray-curtis")


def _distance_array(distance_matrix: pd.DataFrame) -> np.ndarray:
    distances = distance_matrix.to_numpy(dtype=float)
    if distances.shape[0] != distances.shape[1]:
        raise MicrobiomeSuiteError("Distance matrix must be square.")
    if not np.allclose(distances, distances.T):
        raise MicrobiomeSuiteError("Distance matrix must be symmetric.")
    return distances


def _lower_triangle(matrix: np.ndarray) -> np.ndarray:
    return matrix[np.tril_indices(matrix.shape[0], k=-1)]


def _vector_correlation(x: np.ndarray, y: np.ndarray, *, method: str) -> float:
    method = method.lower()
    if method == "spearman":
        x = pd.Series(x).rank(method="average").to_numpy()
        y = pd.Series(y).rank(method="average").to_numpy()
    elif method != "pearson":
        raise MicrobiomeSuiteError("--method must be pearson or spearman.")
    if np.std(x) == 0 or np.std(y) == 0:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])
