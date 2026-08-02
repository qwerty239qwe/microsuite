from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from microsuite.methods._sparcc import estimate_sparcc

FIXTURE_DIR = Path(__file__).resolve().parents[1] / "fixtures" / "sparcc"
RUNNER = Path(__file__).with_name("run_spieceasi_reference.R")
REFERENCE_SEEDS = (10010, 10011, 10012)
ALLOWED_MAE = 0.05300942276035425

pytestmark = pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_SPARCC_REFERENCE") != "1",
    reason="set MICROSUITE_RUN_SPARCC_REFERENCE=1 to run the live SpiecEasi check",
)


def _vgam_available(rscript: str) -> bool:
    probe = subprocess.run(
        [rscript, "-e", 'quit(status = !requireNamespace("VGAM", quietly = TRUE))'],
        capture_output=True,
        text=True,
        check=False,
    )
    return probe.returncode == 0


def _read_matrix(path: Path, expected_features: list[str]) -> np.ndarray:
    frame = pd.read_csv(path, sep="\t", index_col=0)
    assert frame.index.tolist() == expected_features
    assert frame.columns.tolist() == expected_features
    values = frame.to_numpy(dtype=float)
    assert values.shape == (len(expected_features), len(expected_features))
    assert np.isfinite(values).all()
    np.testing.assert_allclose(values, values.T, rtol=0.0, atol=1e-14)
    np.testing.assert_allclose(np.diag(values), 1.0, rtol=0.0, atol=1e-14)
    return values


def _off_diagonal_mae(left: np.ndarray, right: np.ndarray) -> float:
    upper = np.triu_indices_from(left, k=1)
    return float(np.mean(np.abs(left[upper] - right[upper])))


def test_live_spieceasi_matches_microsuite_and_frozen_evidence(tmp_path: Path) -> None:
    rscript = shutil.which("Rscript")
    if rscript is None:
        pytest.skip("Rscript is required for the live SpiecEasi check")
    if not _vgam_available(rscript):
        pytest.skip("Rscript with the VGAM package is required")

    source_dir_value = os.environ.get("SPIECEASI_SOURCE_DIR")
    if not source_dir_value:
        pytest.skip("set SPIECEASI_SOURCE_DIR to the pinned SpiecEasi checkout")
    source_dir = Path(source_dir_value)
    assert source_dir.is_dir(), f"SPIECEASI_SOURCE_DIR is not a directory: {source_dir}"

    subprocess.run(
        [
            rscript,
            str(RUNNER),
            str(source_dir),
            str(FIXTURE_DIR),
            str(tmp_path),
        ],
        check=True,
    )

    truth_frame = pd.read_csv(FIXTURE_DIR / "latent_correlation.tsv", sep="\t", index_col=0)
    features = truth_frame.columns.tolist()
    assert truth_frame.index.tolist() == features
    truth = truth_frame.to_numpy(dtype=float)
    strong_pairs = np.argwhere(np.triu(np.abs(truth) >= 0.7, k=1))
    assert strong_pairs.shape == (3, 2)

    for dataset in ("dense", "zero"):
        counts_frame = pd.read_csv(FIXTURE_DIR / f"{dataset}_counts.tsv", sep="\t", index_col=0)
        assert counts_frame.columns.tolist() == features
        candidate = estimate_sparcc(counts_frame.to_numpy(dtype=float), seed=10010).correlation
        live_reference = _read_matrix(tmp_path / f"{dataset}_reference_cor.tsv", features)
        frozen_references = np.stack(
            [
                _read_matrix(
                    FIXTURE_DIR / f"{dataset}_reference_cor_seed_{seed}.tsv",
                    features,
                )
                for seed in REFERENCE_SEEDS
            ]
        )
        frozen_median = np.median(frozen_references, axis=0)

        candidate_mae = _off_diagonal_mae(candidate, live_reference)
        capture_mae = _off_diagonal_mae(live_reference, frozen_median)
        assert candidate_mae <= ALLOWED_MAE, (
            f"{dataset} microsuite-to-live SpiecEasi MAE {candidate_mae} exceeds "
            f"the frozen limit {ALLOWED_MAE}"
        )
        assert capture_mae <= ALLOWED_MAE, (
            f"{dataset} live-to-frozen SpiecEasi MAE {capture_mae} exceeds "
            f"the frozen limit {ALLOWED_MAE}"
        )

        for left, right in strong_pairs:
            expected_sign = np.sign(truth[left, right])
            assert np.sign(candidate[left, right]) == expected_sign
            assert np.sign(live_reference[left, right]) == expected_sign
