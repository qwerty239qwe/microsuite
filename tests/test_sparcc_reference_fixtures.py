from __future__ import annotations

import hashlib
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "sparcc"
COUNT_FEATURES = [f"feature_{index:02d}" for index in range(1, 11)]
INNER_FEATURES = [f"inner_feature_{index:02d}" for index in range(1, 7)]
SEEDS = (10010, 10011, 10012)
ALLOWED_MAE = 0.05300942276035425

# These are evidence-integrity checks, not candidate-implementation
# expectations. README.md is intentionally excluded because it records them.
EXPECTED_SHA256 = {
    "generate_inputs.py": "1f807d24b6d5279097fa73de9c6ed19dc67dbb044b5649843aa347940bb4cd63",
    "capture_reference.R": "aea4b40cde4d089bfd4b23a6a8704a95811dfdc7ab9acd28c59405d734dbf9a3",
    "dense_counts.tsv": "6f3f5d9591ab6a6ad94c7a66948be2dd3cf1ee0e52862a3b2db85b57707d76c8",
    "zero_counts.tsv": "741273ed55fd75e5f07bcce6ff7b3d3efee4342faecd828e1a281b5974c940da",
    "inner_compositions.tsv": "dd0e611b85863fe88a7f2e66a9c30b014032d9e9c530bb1b97ffe5f4ecbd79d9",
    "latent_correlation.tsv": "b1aa6e2b318cb94efa622a7328492626241390b15553651ef73cd5ac9101be38",
    "inner_initial_reference_cor.tsv": (
        "ec12047a7c9954f3be4e421714158120143faf86c872805ca09dbdb56d0703b9"
    ),
    "inner_reference_cor.tsv": "bac0fb0f1e291a7fd0418598101717ee8683ebda1a30a4f2b1316af06ba2750f",
    "dense_reference_cor_seed_10010.tsv": (
        "9f4cf0eaf3715bb3d9620ba1973754c880aefbcaf3ccdfb7449534e750dda007"
    ),
    "dense_reference_cor_seed_10011.tsv": (
        "2e84673c572f48de0bf02feeb1d1984abed3a138481acb51c335104d9a2f1e42"
    ),
    "dense_reference_cor_seed_10012.tsv": (
        "7edf4d3ed2c0ab10f613e4a59410509a4852b503b1b855073a9aaeece4e8a0fb"
    ),
    "zero_reference_cor_seed_10010.tsv": (
        "0661bc47815604fb91edd7945801bde23974dedcb12c6e559804254510f9c335"
    ),
    "zero_reference_cor_seed_10011.tsv": (
        "ec46876418136c32487ac7d4dab2ff4e7286e8645563c83c5987c4213bc4dbd6"
    ),
    "zero_reference_cor_seed_10012.tsv": (
        "2fe541ab65bfab9aa81875c344b3db4a57df729428a2dbbd66f979a5a254af7b"
    ),
}


def _read_table(filename: str) -> pd.DataFrame:
    return pd.read_csv(FIXTURE_DIR / filename, sep="\t", index_col=0)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _off_diagonal_mae(left: np.ndarray, right: np.ndarray) -> float:
    upper = np.triu_indices_from(left, k=1)
    return float(np.mean(np.abs(left[upper] - right[upper])))


def test_input_tables_have_frozen_shapes_order_and_domains() -> None:
    dense = _read_table("dense_counts.tsv")
    zero = _read_table("zero_counts.tsv")
    inner = _read_table("inner_compositions.tsv")
    truth = _read_table("latent_correlation.tsv")

    for counts in (dense, zero):
        assert counts.shape == (400, 10)
        assert counts.columns.tolist() == COUNT_FEATURES
        assert counts.index.tolist() == [f"sample_{index:03d}" for index in range(1, 401)]
        values = counts.to_numpy(dtype=float)
        assert np.isfinite(values).all()
        assert (values >= 0).all()
        assert np.equal(values, np.floor(values)).all()
        assert (values.sum(axis=0) > 0).all()
        assert (values.sum(axis=1) > 0).all()

    assert 0.18 <= np.mean(zero.to_numpy() == 0) <= 0.20
    assert inner.shape == (24, 6)
    assert inner.columns.tolist() == INNER_FEATURES
    assert inner.index.tolist() == [f"inner_sample_{index:02d}" for index in range(1, 25)]
    assert np.isfinite(inner.to_numpy()).all()
    assert (inner.to_numpy() > 0).all()
    np.testing.assert_allclose(inner.sum(axis=1), 1.0, rtol=0.0, atol=2e-15)

    assert truth.shape == (10, 10)
    assert truth.index.tolist() == COUNT_FEATURES
    assert truth.columns.tolist() == COUNT_FEATURES
    np.testing.assert_allclose(truth, truth.T, rtol=0.0, atol=1e-15)
    np.testing.assert_allclose(np.diag(truth), 1.0, rtol=0.0, atol=1e-15)
    assert np.linalg.eigvalsh(truth.to_numpy()).min() > 0.0


@pytest.mark.parametrize("dataset", ["dense", "zero"])
@pytest.mark.parametrize("seed", SEEDS)
def test_outer_reference_matrices_are_valid(dataset: str, seed: int) -> None:
    matrix = _read_table(f"{dataset}_reference_cor_seed_{seed}.tsv")
    assert matrix.shape == (10, 10)
    assert matrix.index.tolist() == COUNT_FEATURES
    assert matrix.columns.tolist() == COUNT_FEATURES
    values = matrix.to_numpy(dtype=float)
    assert np.isfinite(values).all()
    np.testing.assert_allclose(values, values.T, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(np.diag(values), 1.0, rtol=0.0, atol=1e-14)
    assert np.max(np.abs(values)) <= 1.0


@pytest.mark.parametrize(
    "filename",
    ["inner_initial_reference_cor.tsv", "inner_reference_cor.tsv"],
)
def test_inner_reference_matrices_are_valid(filename: str) -> None:
    matrix = _read_table(filename)
    assert matrix.shape == (6, 6)
    assert matrix.index.tolist() == INNER_FEATURES
    assert matrix.columns.tolist() == INNER_FEATURES
    values = matrix.to_numpy(dtype=float)
    assert np.isfinite(values).all()
    np.testing.assert_allclose(values, values.T, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(np.diag(values), 1.0, rtol=0.0, atol=1e-14)
    assert np.max(np.abs(values)) <= 1.0


def test_acceptance_threshold_is_derived_only_from_reference_variability() -> None:
    maximum_reference_seed_mae = 0.0
    for dataset in ("dense", "zero"):
        references = {
            seed: _read_table(f"{dataset}_reference_cor_seed_{seed}.tsv").to_numpy()
            for seed in SEEDS
        }
        for left_index, left_seed in enumerate(SEEDS):
            for right_seed in SEEDS[left_index + 1 :]:
                maximum_reference_seed_mae = max(
                    maximum_reference_seed_mae,
                    _off_diagonal_mae(references[left_seed], references[right_seed]),
                )

    derived_allowed_mae = max(0.02, 5 * maximum_reference_seed_mae)
    assert maximum_reference_seed_mae == pytest.approx(0.01060188455207085, abs=1e-16)
    assert derived_allowed_mae == pytest.approx(ALLOWED_MAE, abs=1e-16)

    readme = (FIXTURE_DIR / "README.md").read_text(encoding="utf-8")
    recorded = re.search(r"allowed_mae = max.*?\n\s*= ([0-9.]+)", readme)
    assert recorded is not None
    assert float(recorded.group(1)) == pytest.approx(derived_allowed_mae, abs=1e-16)


def test_provenance_and_hashes_are_frozen() -> None:
    readme = (FIXTURE_DIR / "README.md").read_text(encoding="utf-8")
    assert "generated evidence; never hand-edit" in readme.lower()
    assert "faed6a4476fe0a8dc701ea15cbdfe98d56ce6704" in readme
    assert "R version 4.6.0 (2026-04-24)" in readme
    assert "VGAM" in readme and "1.1-14" in readme
    assert "samples by features" in readme

    assert EXPECTED_SHA256
    for filename, expected_hash in EXPECTED_SHA256.items():
        assert _sha256(FIXTURE_DIR / filename) == expected_hash
        assert f"`{filename}`" in readme
        assert f"`{expected_hash}`" in readme


def test_python_inputs_regenerate_byte_for_byte(tmp_path: Path) -> None:
    subprocess.run(
        [
            sys.executable,
            str(FIXTURE_DIR / "generate_inputs.py"),
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
    )
    for filename in (
        "dense_counts.tsv",
        "zero_counts.tsv",
        "inner_compositions.tsv",
        "latent_correlation.tsv",
    ):
        assert _sha256(tmp_path / filename) == _sha256(FIXTURE_DIR / filename)
