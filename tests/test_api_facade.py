from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import microsuite.api as api
from microsuite._errors import MicrobiomeSuiteError
from microsuite.api import (
    abundance_table,
    alpha_diversity,
    beta_diversity,
    normalize_table,
    pcoa,
    qc,
    rarefy_table,
    read_table,
    shared_taxa_table,
    write_table,
)
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def test_python_sdk_facade_table_roundtrip_and_ecology(tmp_path: Path) -> None:
    table = tmp_path / "table.h5ad"
    write_table(fixture_adata(), table)
    adata = read_table(table)

    relative = normalize_table(adata, method="relative")
    abundance = abundance_table(adata, level="genus")
    shared = shared_taxa_table(adata, level="genus", group="body_site")
    rarefied = rarefy_table(adata, depth=10, seed=1)
    alpha = alpha_diversity(adata, metric="shannon")
    beta = beta_diversity(adata, metric="bray-curtis")
    coords = pcoa(beta, dimensions=2)

    assert np.allclose(np.asarray(relative.X).sum(axis=1), 1.0)
    assert isinstance(abundance, pd.DataFrame)
    assert "Lactobacillus" in abundance.columns
    assert shared["taxon"].tolist()
    assert np.asarray(rarefied.X).sum(axis=1).tolist() == [10.0, 10.0, 10.0, 10.0]
    assert "shannon" in alpha.columns
    assert beta.shape == (4, 4)
    assert coords.columns.tolist() == [
        "sample_id",
        "PC1",
        "PC2",
        "PC1_variance",
        "PC2_variance",
    ]


def test_python_sdk_facade_public_exports_include_qc() -> None:
    assert "qc" in api.__all__
    assert "functional_profile" in api.__all__
    assert "beta_significance" in api.__all__
    assert "mantel_test" in api.__all__
    assert "gamma_diversity" in api.__all__
    assert "constrained_ordination" in api.__all__
    assert "network" in api.__all__
    assert "ml_classify" in api.__all__
    assert "longitudinal" in api.__all__


def test_python_sdk_facade_exposes_fastqc(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    read = tmp_path / "sample_R1.fastq.gz"
    read.touch()
    calls: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "fastqc" if name == "fastqc" else None)

    def fake_run(command: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        calls.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr("subprocess.run", fake_run)

    qc(
        backend="fastqc",
        inputs=[read],
        output_dir=tmp_path / "qc",
        threads=4,
        extract=True,
    )

    assert calls == [
        [
            "fastqc",
            "--outdir",
            str(tmp_path / "qc"),
            "--threads",
            "4",
            "--extract",
            str(read),
        ]
    ]


def test_python_sdk_facade_fastqc_reports_missing_binary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    read = tmp_path / "sample_R1.fastq.gz"
    read.touch()
    monkeypatch.setattr("shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="FastQC requires"):
        qc(backend="fastqc", inputs=[read], output_dir=tmp_path / "qc")
