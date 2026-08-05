from __future__ import annotations

import shutil
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity.vegan import _invoke_vegan, vegan_beta_significance

FIXTURE = Path(__file__).parent.parent / "containers" / "r-ecology" / "smoke"


def fixture_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    distance = pd.read_csv(FIXTURE / "distance.tsv", sep="\t", index_col=0)
    metadata = pd.read_csv(FIXTURE / "metadata.tsv", sep="\t", index_col=0)
    return distance, metadata


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_vegan_adonis2_and_anosim2_outputs_are_tidy() -> None:
    distance, metadata = fixture_inputs()

    adonis = vegan_beta_significance(distance, metadata, formula="group", permutations=19)

    assert adonis.loc[0, "method"] == "adonis2"
    assert {"term", "r_squared", "f_value", "p_value"}.issubset(adonis.columns)

    anosim = vegan_beta_significance(
        distance,
        metadata,
        column="group",
        method="anosim2",
        permutations=19,
        seed=4,
    )
    assert anosim.loc[0, "method"] == "anosim2"
    assert anosim.loc[0, "n_samples"] == 6
    assert 0 <= anosim.loc[0, "p_value"] <= 1


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_vegan_formula_and_strata_are_reproducible() -> None:
    distance, metadata = fixture_inputs()

    first = vegan_beta_significance(
        distance,
        metadata,
        formula="group",
        strata="block",
        method="adonis2",
        permutations=19,
        seed=7,
    )
    second = vegan_beta_significance(
        distance,
        metadata,
        formula="group",
        strata="block",
        method="adonis2",
        permutations=19,
        seed=7,
    )

    pd.testing.assert_frame_equal(first, second)
    assert first.loc[0, "permutation_scheme"] == "blocked"


@pytest.mark.skipif(shutil.which("Rscript") is None, reason="Rscript is not installed")
def test_vegan_nested_formula_strata_and_permutations_work_together() -> None:
    distance, metadata = fixture_inputs()

    result = vegan_beta_significance(
        distance,
        metadata,
        formula="group / block",
        strata="group:block",
        method="adonis2",
        permutations=19,
        seed=11,
    )

    assert result.loc[0, "formula"] == "group / block"
    assert result.loc[0, "strata"] == "group:block"
    assert result["permutations"].eq(19).all()
    assert result["permutation_scheme"].eq("blocked").all()


def test_vegan_validates_formula_and_anosim_shape() -> None:
    distance, metadata = fixture_inputs()

    with pytest.raises(MicrobiomeSuiteError, match="not found"):
        vegan_beta_significance(distance, metadata, formula="missing_group")
    with pytest.raises(MicrobiomeSuiteError, match="one grouping column"):
        vegan_beta_significance(distance, metadata, formula="group + block", method="anosim2")
    with pytest.raises(MicrobiomeSuiteError, match="does not support lme4-style"):
        vegan_beta_significance(distance, metadata, formula="group + (1 | group:block)")


def test_vegan_local_runtime_reports_missing_rscript(monkeypatch) -> None:
    distance, metadata = fixture_inputs()
    monkeypatch.setattr("microsuite.diversity.vegan.shutil.which", lambda name: None)

    with pytest.raises(MicrobiomeSuiteError, match="Rscript"):
        vegan_beta_significance(distance, metadata, formula="group")


def test_vegan_docker_command_uses_ecology_image_and_uid(tmp_path: Path, monkeypatch) -> None:
    distance = tmp_path / "distance.tsv"
    metadata = tmp_path / "metadata.tsv"
    output = tmp_path / "result.tsv"
    distance.write_text("\tS0\nS0\t0\n", encoding="utf-8")
    metadata.write_text("\tgroup\nS0\tA\n", encoding="utf-8")
    commands: list[list[str]] = []

    monkeypatch.setattr("shutil.which", lambda name: "/usr/bin/docker")
    monkeypatch.setattr("microsuite.diversity.vegan.host_user_spec", lambda: "1000:1000")

    def fake_run(command: list[str], **kwargs: object) -> SimpleNamespace:
        commands.append(command)
        output.write_text("backend\tmethod\nvegan\tadonis2\n", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("microsuite.diversity.vegan.run_command", fake_run)
    monkeypatch.setattr(
        "microsuite.diversity.vegan.resolve_image_digest", lambda *args: "sha256:test"
    )

    _invoke_vegan(
        method="adonis2",
        distance_path=distance,
        metadata_path=metadata,
        formula="group",
        strata=None,
        permutations=9,
        seed=0,
        output_path=output,
        runtime="docker",
        image="example/r-ecology:test",
        engine="docker",
        run_dir=None,
        timeout=None,
        sidecar_dir=tmp_path,
    )

    assert commands[0][:4] == ["docker", "run", "--rm", "--user"]
    assert "example/r-ecology:test" in commands[0]
    assert (tmp_path / "vegan_container.json").exists()
