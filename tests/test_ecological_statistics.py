from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from typer.testing import CliRunner

from microsuite.api import (
    beta_significance,
    beta_turnover,
    constrained_ordination,
    gamma_diversity,
    mantel_test,
    taxa_turnover,
)
from microsuite.cli.app import app
from microsuite.diversity.beta import beta_diversity
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_adata():
    return read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")


def test_beta_significance_permanova_permdisp_and_anosim() -> None:
    adata = fixture_adata()
    beta = beta_diversity(adata, "bray-curtis")
    metadata = pd.DataFrame(adata.obs)

    permanova = beta_significance(
        beta,
        metadata,
        column="body_site",
        method="permanova",
        permutations=19,
        seed=1,
    )
    permdisp = beta_significance(
        beta,
        metadata,
        column="body_site",
        method="permdisp",
        permutations=19,
        seed=1,
    )
    anosim = beta_significance(
        beta,
        metadata,
        column="body_site",
        method="anosim",
        permutations=19,
        seed=1,
    )

    assert permanova.loc[0, "method"] == "permanova"
    assert permanova.loc[0, "n_groups"] == 2
    assert 0 <= permanova.loc[0, "p_value"] <= 1
    assert permdisp.loc[0, "method"] == "permdisp"
    assert permdisp.loc[0, "n_groups"] == 2
    assert permdisp.loc[0, "f_value"] >= 0
    assert 0 <= permdisp.loc[0, "p_value"] <= 1
    assert permdisp.loc[0, "permutation_scheme"] == "unrestricted"
    assert anosim.loc[0, "method"] == "anosim"
    assert "r" in anosim.columns


def test_mantel_test_correlates_distance_matrices() -> None:
    adata = fixture_adata()
    bray = beta_diversity(adata, "bray-curtis")
    jaccard = beta_diversity(adata, "jaccard")

    result = mantel_test(bray, jaccard, method="spearman", permutations=19, seed=2)

    assert result.loc[0, "method"] == "mantel-spearman"
    assert result.loc[0, "n_samples"] == adata.n_obs
    assert -1 <= result.loc[0, "r"] <= 1


def test_gamma_and_turnover_outputs() -> None:
    adata = fixture_adata()

    gamma = gamma_diversity(adata, group="body_site", metric="observed_features")
    beta = beta_turnover(adata, level="genus")
    taxa = taxa_turnover(adata, group="body_site", level="genus")

    assert gamma["group"].tolist() == ["gut", "tongue"]
    assert gamma["observed_features"].min() > 0
    assert {"sample_a", "sample_b", "turnover_component"}.issubset(beta.columns)
    assert taxa.loc[0, "group_a"] == "gut"
    assert taxa.loc[0, "group_b"] == "tongue"
    assert 0 <= taxa.loc[0, "jaccard_turnover"] <= 1


def test_constrained_ordination_methods_return_coordinates() -> None:
    adata = fixture_adata()

    for method in ["rda", "db-rda", "cca"]:
        result = constrained_ordination(
            adata,
            constraints=["body_site"],
            method=method,
            dimensions=2,
        )

        assert result["sample_id"].tolist() == adata.obs_names.astype(str).tolist()
        assert result["method"].unique().tolist() == [method]
        assert np.isfinite(result[["Axis1", "Axis2"]].to_numpy()).all()


def test_diversity_cli_beta_significance_and_gamma(tmp_path: Path) -> None:
    adata = fixture_adata()
    table = tmp_path / "table.h5ad"
    distance = tmp_path / "beta.tsv"
    metadata = tmp_path / "metadata.tsv"
    beta_output = tmp_path / "permanova.tsv"
    gamma_output = tmp_path / "gamma.tsv"

    write_h5ad(adata, table)
    beta_diversity(adata, "bray-curtis").to_csv(distance, sep="\t")
    pd.DataFrame(adata.obs).to_csv(metadata, sep="\t")

    runner = CliRunner()
    beta_result = runner.invoke(
        app,
        [
            "diversity",
            "beta-significance",
            str(distance),
            "--metadata",
            str(metadata),
            "--column",
            "body_site",
            "--permutations",
            "9",
            "-o",
            str(beta_output),
        ],
    )
    gamma_result = runner.invoke(
        app,
        [
            "diversity",
            "gamma",
            str(table),
            "--group",
            "body_site",
            "-o",
            str(gamma_output),
        ],
    )

    assert beta_result.exit_code == 0, beta_result.stdout
    assert gamma_result.exit_code == 0, gamma_result.stdout
    assert pd.read_csv(beta_output, sep="\t").loc[0, "method"] == "permanova"
    assert pd.read_csv(gamma_output, sep="\t").shape[0] == 2
