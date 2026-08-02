#!/usr/bin/env python3
"""Generate deterministic, implementation-independent SparCC fixture inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

SEED = 10010
N_SAMPLES = 400
COUNT_FEATURES = tuple(f"feature_{index:02d}" for index in range(1, 11))
INNER_FEATURES = tuple(f"inner_feature_{index:02d}" for index in range(1, 7))


def latent_correlation() -> np.ndarray:
    """Return a sparse positive-definite latent log-abundance correlation."""
    loadings = np.array(
        [
            [0.90, 0.00, 0.00],
            [0.85, 0.00, 0.00],
            [-0.85, 0.00, 0.00],
            [-0.45, 0.00, 0.00],
            [0.00, 0.75, 0.00],
            [0.00, 0.55, 0.00],
            [0.00, -0.65, 0.00],
            [0.00, 0.00, 0.70],
            [0.00, 0.00, -0.55],
            [0.00, 0.00, 0.45],
        ],
        dtype=np.float64,
    )
    unique_variance = 1.0 - np.einsum("ij,ij->i", loadings, loadings)
    correlation = loadings @ loadings.T + np.diag(unique_variance)
    if not np.allclose(np.diag(correlation), 1.0):
        raise RuntimeError("latent correlation does not have a unit diagonal")
    if np.linalg.eigvalsh(correlation).min() <= 0.0:
        raise RuntimeError("latent correlation is not positive definite")
    return correlation


def _softmax_rows(values: np.ndarray) -> np.ndarray:
    shifted = values - values.max(axis=1, keepdims=True)
    weights = np.exp(shifted)
    return weights / weights.sum(axis=1, keepdims=True)


def _draw_counts(
    rng: np.random.Generator,
    correlation: np.ndarray,
    depths: np.ndarray,
) -> np.ndarray:
    means = np.linspace(-0.30, 0.30, len(COUNT_FEATURES))
    log_abundances = rng.multivariate_normal(
        mean=means,
        cov=correlation,
        size=N_SAMPLES,
        check_valid="raise",
    )
    probabilities = _softmax_rows(log_abundances)
    return np.vstack(
        [
            rng.multinomial(int(depth), probability)
            for depth, probability in zip(depths, probabilities, strict=True)
        ]
    )


def _write_table(
    path: Path,
    values: np.ndarray,
    row_names: tuple[str, ...],
    column_names: tuple[str, ...],
    *,
    integer: bool = False,
) -> None:
    if values.shape != (len(row_names), len(column_names)):
        raise ValueError(f"shape mismatch while writing {path}")
    formatter = (
        (lambda value: str(int(value))) if integer else (lambda value: format(float(value), ".17g"))
    )
    lines = ["\t".join(("id", *column_names))]
    for row_name, row in zip(row_names, values, strict=True):
        lines.append("\t".join((row_name, *(formatter(value) for value in row))))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    correlation = latent_correlation()

    dense_depths = np.rint(
        np.exp(rng.uniform(np.log(500.0), np.log(20_000.0), size=N_SAMPLES))
    ).astype(np.int64)
    dense_counts = _draw_counts(rng, correlation, dense_depths)

    zero_depths = rng.integers(250, 1_001, size=N_SAMPLES)
    zero_counts = _draw_counts(rng, correlation, zero_depths)
    target_zero_cells = int(round(0.19 * zero_counts.size))
    positive_indices = np.flatnonzero(zero_counts)
    additional_zeros = target_zero_cells - int(np.count_nonzero(zero_counts == 0))
    tie_breakers = rng.random(positive_indices.size)
    detection_order = np.lexsort((tie_breakers, zero_counts.flat[positive_indices]))
    dropout_indices = positive_indices[detection_order[:additional_zeros]]
    zero_counts.flat[dropout_indices] = 0

    if np.any(dense_counts.sum(axis=0) == 0) or np.any(dense_counts.sum(axis=1) == 0):
        raise RuntimeError("dense fixture contains an all-zero row or feature")
    if np.any(zero_counts.sum(axis=0) == 0) or np.any(zero_counts.sum(axis=1) == 0):
        raise RuntimeError("zero fixture contains an all-zero row or feature")

    inner_compositions = rng.dirichlet(np.linspace(1.5, 4.0, len(INNER_FEATURES)), size=24)

    sample_names = tuple(f"sample_{index:03d}" for index in range(1, N_SAMPLES + 1))
    inner_sample_names = tuple(f"inner_sample_{index:02d}" for index in range(1, 25))
    _write_table(
        output_dir / "dense_counts.tsv",
        dense_counts,
        sample_names,
        COUNT_FEATURES,
        integer=True,
    )
    _write_table(
        output_dir / "zero_counts.tsv",
        zero_counts,
        sample_names,
        COUNT_FEATURES,
        integer=True,
    )
    _write_table(
        output_dir / "inner_compositions.tsv",
        inner_compositions,
        inner_sample_names,
        INNER_FEATURES,
    )
    _write_table(
        output_dir / "latent_correlation.tsv",
        correlation,
        COUNT_FEATURES,
        COUNT_FEATURES,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="fixture directory (default: directory containing this script)",
    )
    args = parser.parse_args()
    generate(args.output_dir.resolve())


if __name__ == "__main__":
    main()
