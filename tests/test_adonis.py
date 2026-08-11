from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from typer.testing import CliRunner

from microsuite._errors import MicrobiomeSuiteError
from microsuite.api import adonis2
from microsuite.cli.app import app
from microsuite.diversity.adonis import (
    PermutationScheme,
    _achievable_plot_permutations,
    parse_formula,
    permutation_indices,
    scheme_from_formula,
)
from microsuite.diversity.beta import beta_diversity
from microsuite.diversity.vegan import vegan_beta_significance
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def fixture_beta_and_metadata() -> tuple[pd.DataFrame, pd.DataFrame]:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    return beta_diversity(adata, "bray-curtis"), pd.DataFrame(adata.obs)


def test_parse_formula_expands_terms_and_random_effects() -> None:
    formula = parse_formula(
        "counts ~ disease_status + accession + accession:timepoint + (1 | accession:subject_id)"
    )
    assert [term.label for term in formula.terms] == [
        "disease_status",
        "accession",
        "accession:timepoint",
    ]
    assert formula.random_groups == (("accession", "subject_id"),)
    assert formula.columns == ("disease_status", "accession", "timepoint", "subject_id")


def test_parse_formula_expands_crossing() -> None:
    assert [term.label for term in parse_formula("~ a*b").terms] == ["a", "b", "a:b"]


@pytest.mark.parametrize(
    "formula",
    ["counts ~", "~ (a | b)", "~ (1 | )", "~ a + (a"],
)
def test_parse_formula_rejects_malformed_input(formula: str) -> None:
    with pytest.raises(MicrobiomeSuiteError):
        parse_formula(formula)


def test_adonis2_decomposition_is_additive() -> None:
    beta, metadata = fixture_beta_and_metadata()
    result = adonis2(beta, metadata, formula="d ~ body_site + subject", permutations=0, seed=1)
    terms = result[~result["term"].isin({"Residual", "Total"})]
    total = float(result.loc[result["term"] == "Total", "sum_of_squares"].iloc[0])
    residual = float(result.loc[result["term"] == "Residual", "sum_of_squares"].iloc[0])

    assert terms["sum_of_squares"].sum() + residual == pytest.approx(total)
    assert result["r2"].iloc[:-1].sum() == pytest.approx(1.0)
    n_samples = int(result["n_samples"].iloc[0])
    assert int(result.loc[result["term"] == "Total", "df"].iloc[0]) == n_samples - 1
    assert (
        terms["df"].sum() + int(result.loc[result["term"] == "Residual", "df"].iloc[0])
        == n_samples - 1
    )


def test_adonis2_single_term_matches_one_way_permanova() -> None:
    """With one factor the sequential decomposition reduces to the classic test."""
    from microsuite.diversity.ecology import _permanova_statistic

    beta, metadata = fixture_beta_and_metadata()
    result = adonis2(beta, metadata, formula="d ~ body_site", permutations=0, seed=1)
    samples = [s for s in beta.index.astype(str) if s in set(metadata.index.astype(str))]
    labels = metadata.loc[samples, "body_site"].astype(str).to_numpy()
    expected = _permanova_statistic(beta.loc[samples, samples].to_numpy(dtype=float), labels)

    assert float(result.loc[0, "pseudo_f"]) == pytest.approx(expected, rel=1e-9)


def test_adonis2_permutation_p_values_are_bounded() -> None:
    beta, metadata = fixture_beta_and_metadata()
    result = adonis2(beta, metadata, formula="d ~ body_site", permutations=49, seed=3)
    p_value = float(result.loc[0, "p_value"])
    assert 1 / 50 <= p_value <= 1.0
    assert result.loc[0, "permutation_scheme"] == "free"


def test_adonis2_is_reproducible_under_a_seed() -> None:
    beta, metadata = fixture_beta_and_metadata()
    kwargs = dict(formula="d ~ body_site", permutations=49, seed=11)
    first = adonis2(beta, metadata, **kwargs)
    second = adonis2(beta, metadata, **kwargs)
    pd.testing.assert_frame_equal(first, second)


def test_adonis2_random_term_restricts_permutation() -> None:
    beta, metadata = fixture_beta_and_metadata()
    result = adonis2(
        beta,
        metadata,
        formula="d ~ body_site + (1 | subject)",
        permutations=49,
        seed=5,
    )
    assert result.loc[0, "permutation_scheme"] == "plots-free"


def test_adonis2_rejects_a_fully_aliased_term() -> None:
    beta, metadata = fixture_beta_and_metadata()
    metadata = metadata.copy()
    metadata["copy_of_body_site"] = metadata["body_site"]
    with pytest.raises(MicrobiomeSuiteError, match="aliased"):
        adonis2(
            beta,
            metadata,
            formula="d ~ body_site + copy_of_body_site",
            permutations=0,
        )


def test_adonis2_reports_missing_columns() -> None:
    beta, metadata = fixture_beta_and_metadata()
    with pytest.raises(MicrobiomeSuiteError, match="not found"):
        adonis2(beta, metadata, formula="d ~ nope", permutations=0)


def test_permutation_indices_respect_blocks_and_plots() -> None:
    rng = np.random.default_rng(0)
    blocks = np.repeat([0, 1], 12)
    plots = np.repeat(np.arange(8), 3)
    scheme = PermutationScheme(blocks=blocks, plots=plots, within="free")

    for _ in range(25):
        perm = permutation_indices(scheme, 24, rng)
        assert np.array_equal(np.sort(perm), np.arange(24))
        assert np.array_equal(blocks[perm], blocks)
        # Each plot is filled from exactly one donor plot.
        donors = pd.DataFrame({"target": plots, "donor": plots[perm]})
        assert (donors.groupby("target")["donor"].nunique() == 1).all()


def test_permutation_indices_only_swap_equal_sized_plots() -> None:
    rng = np.random.default_rng(1)
    plots = np.array([0, 0, 0, 1, 1, 2])
    scheme = PermutationScheme(plots=plots, within="none")
    sizes = pd.Series(plots).value_counts()

    for _ in range(25):
        perm = permutation_indices(scheme, 6, rng)
        assert np.array_equal(np.sort(perm), np.arange(6))
        donors = pd.DataFrame({"target": plots, "donor": plots[perm]})
        for target, donor in donors.groupby("target")["donor"].first().items():
            assert sizes[target] == sizes[donor]


def test_permutation_indices_within_none_keeps_sample_order() -> None:
    rng = np.random.default_rng(2)
    plots = np.repeat(np.arange(4), 2)
    scheme = PermutationScheme(plots=plots, within="none")
    perm = permutation_indices(scheme, 8, rng)
    # Within a plot the donor samples keep their relative order.
    for plot in range(4):
        positions = np.flatnonzero(plots == plot)
        assert np.all(np.diff(perm[positions]) == 1)


def test_scheme_from_formula_nests_subjects_inside_the_outer_factor() -> None:
    metadata = pd.DataFrame(
        {
            "accession": ["A", "A", "B", "B"],
            "subject_id": ["1", "1", "1", "1"],
        }
    )
    formula = parse_formula("~ accession + (1 | accession:subject_id)")
    scheme = scheme_from_formula(formula, metadata)
    # subject_id repeats across studies, so the key must include accession.
    assert len(set(scheme.plots.tolist())) == 2
    assert scheme.blocks is not None


def test_scheme_from_formula_rejects_multiple_random_terms() -> None:
    formula = parse_formula("~ a + (1 | b) + (1 | c)")
    with pytest.raises(MicrobiomeSuiteError, match="Only one random-intercept"):
        scheme_from_formula(formula, pd.DataFrame({"a": ["x"], "b": ["y"], "c": ["z"]}))


def test_adonis_cli_writes_a_table(tmp_path: Path) -> None:
    beta, metadata = fixture_beta_and_metadata()
    matrix_path = tmp_path / "dist.tsv"
    metadata_path = tmp_path / "meta.tsv"
    output = tmp_path / "adonis.tsv"
    beta.to_csv(matrix_path, sep="\t")
    metadata.to_csv(metadata_path, sep="\t")

    result = CliRunner().invoke(
        app,
        [
            "diversity",
            "adonis",
            str(matrix_path),
            "--metadata",
            str(metadata_path),
            "--formula",
            "d ~ body_site + (1 | subject)",
            "--permutations",
            "19",
            "--output",
            str(output),
        ],
    )
    assert result.exit_code == 0, result.output
    written = pd.read_csv(output, sep="\t")
    assert list(written["term"])[-2:] == ["Residual", "Total"]
    assert written["permutation_scheme"].iloc[0] == "plots-free"


def test_adonis2_marginal_is_order_invariant() -> None:
    """The point of marginal SS: shared variance goes to nobody, so order stops mattering."""
    beta, metadata = fixture_beta_and_metadata()
    forward = adonis2(
        beta, metadata, formula="d ~ body_site + subject", permutations=0, by="margin"
    )
    reverse = adonis2(
        beta, metadata, formula="d ~ subject + body_site", permutations=0, by="margin"
    )
    for term in ("body_site", "subject"):
        left = forward.loc[forward["term"] == term].iloc[0]
        right = reverse.loc[reverse["term"] == term].iloc[0]
        assert left["df"] == right["df"]
        assert left["sum_of_squares"] == pytest.approx(right["sum_of_squares"])
        assert left["pseudo_f"] == pytest.approx(right["pseudo_f"])


def test_adonis2_rejects_an_unknown_by() -> None:
    beta, metadata = fixture_beta_and_metadata()
    with pytest.raises(MicrobiomeSuiteError, match="terms.*margin"):
        adonis2(beta, metadata, formula="d ~ body_site", permutations=0, by="type-iii")


def confounded_fixture() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Synthetic study/disease design where the two factors overlap.

    Studies A and B are single-disease, C and D are mixed -- the same shape as a
    cross-study meta-analysis, and the shape that makes term order matter.
    """
    rng = np.random.default_rng(0)
    rows = []
    for study, diseases in (
        ("A", ["healthy"] * 8),
        ("B", ["case"] * 8),
        ("C", ["healthy"] * 4 + ["case"] * 4),
        ("D", ["healthy"] * 4 + ["case"] * 4),
    ):
        for index, disease in enumerate(diseases):
            rows.append(
                {
                    "sample_id": f"{study}{index:02d}",
                    "study": study,
                    "disease": disease,
                    "subject": f"{study}-s{index // 2}",
                }
            )
    metadata = pd.DataFrame(rows).set_index("sample_id")

    # Abundances carry a large study shift and a smaller disease shift, so the
    # two factors explain overlapping variance.
    studies = sorted(metadata["study"].unique())
    counts = rng.lognormal(size=(len(metadata), 25))
    for position, study in enumerate(studies):
        rows_for_study = (metadata["study"] == study).to_numpy()
        counts[rows_for_study, position] += 4.0
    counts[(metadata["disease"] == "case").to_numpy(), :3] += 1.5
    proportions = counts / counts.sum(axis=1, keepdims=True)

    difference = np.abs(proportions[:, None, :] - proportions[None, :, :]).sum(axis=2)
    totals = proportions.sum(axis=1)
    beta = pd.DataFrame(
        difference / (totals[:, None] + totals[None, :]),
        index=metadata.index,
        columns=metadata.index,
    )
    return beta, metadata


def test_adonis2_sequential_is_order_dependent_when_terms_are_confounded() -> None:
    beta, metadata = confounded_fixture()
    forward = adonis2(beta, metadata, formula="d ~ disease + study", permutations=0)
    reverse = adonis2(beta, metadata, formula="d ~ study + disease", permutations=0)
    first = float(forward.loc[forward["term"] == "disease", "sum_of_squares"].iloc[0])
    second = float(reverse.loc[reverse["term"] == "disease", "sum_of_squares"].iloc[0])
    # Written first, `disease` also collects the variance it shares with `study`.
    assert first > second

    # The pair explains the same total either way; only the split changes.
    forward_pair = forward[forward["term"].isin({"disease", "study"})]["sum_of_squares"].sum()
    reverse_pair = reverse[reverse["term"].isin({"disease", "study"})]["sum_of_squares"].sum()
    assert forward_pair == pytest.approx(reverse_pair)


def test_adonis2_marginal_recovers_the_unique_contribution() -> None:
    beta, metadata = confounded_fixture()
    reverse = adonis2(beta, metadata, formula="d ~ study + disease", permutations=0)
    marginal = adonis2(beta, metadata, formula="d ~ disease + study", permutations=0, by="margin")
    # For two terms, `disease` entered last already is its marginal contribution.
    assert float(marginal.loc[marginal["term"] == "disease", "sum_of_squares"].iloc[0]) == (
        pytest.approx(float(reverse.loc[reverse["term"] == "disease", "sum_of_squares"].iloc[0]))
    )


def test_adonis2_marginal_reports_a_nested_term_as_zero_df() -> None:
    """`study` carries nothing of its own once `study:subject` is in the model."""
    beta, metadata = confounded_fixture()
    result = adonis2(
        beta, metadata, formula="d ~ study + study:subject", permutations=0, by="margin"
    )
    nested = result.loc[result["term"] == "study"].iloc[0]
    assert int(nested["df"]) == 0
    assert float(nested["sum_of_squares"]) == pytest.approx(0.0)
    assert pd.isna(nested["pseudo_f"])
    assert pd.isna(nested["p_value"])


def test_adonis2_marginal_residual_comes_from_the_full_model() -> None:
    """Marginal SS omit shared variance, so they must not be used to derive the residual."""
    beta, metadata = confounded_fixture()
    sequential = adonis2(beta, metadata, formula="d ~ disease + study", permutations=0)
    marginal = adonis2(beta, metadata, formula="d ~ disease + study", permutations=0, by="margin")

    for column in ("df", "sum_of_squares"):
        assert float(marginal.loc[marginal["term"] == "Residual", column].iloc[0]) == (
            pytest.approx(float(sequential.loc[sequential["term"] == "Residual", column].iloc[0]))
        )

    terms = marginal[~marginal["term"].isin({"Residual", "Total"})]
    residual_ss = float(marginal.loc[marginal["term"] == "Residual", "sum_of_squares"].iloc[0])
    total_ss = float(marginal.loc[marginal["term"] == "Total", "sum_of_squares"].iloc[0])
    # Strictly less, because the confounded portion is credited to neither term.
    assert terms["sum_of_squares"].sum() + residual_ss < total_ss


def test_adonis2_marginal_permutation_p_values_are_bounded() -> None:
    beta, metadata = confounded_fixture()
    result = adonis2(
        beta, metadata, formula="d ~ disease + study", permutations=49, seed=3, by="margin"
    )
    tested = result[~result["term"].isin({"Residual", "Total"})]
    assert ((tested["p_value"] >= 1 / 50) & (tested["p_value"] <= 1.0)).all()


def test_adonis2_marginal_handles_a_zero_df_term_in_the_full_model() -> None:
    """Marginal fitting must tolerate a zero-df term while *building* the model.

    Writing the nested term first -- `study:subject + study` -- leaves `study`
    fully spanned by the time it is reached, so the sequential build underneath
    `build_marginal_design` legitimately produces a zero-width block. That path
    needs `strict=False`; a stricter setting turns a valid formula into an error.
    """
    beta, metadata = confounded_fixture()
    result = adonis2(
        beta, metadata, formula="d ~ study:subject + study", permutations=0, by="margin"
    )
    spanned = result.loc[result["term"] == "study"].iloc[0]
    assert int(spanned["df"]) == 0
    assert float(spanned["sum_of_squares"]) == pytest.approx(0.0)
    # The term that does carry information is unaffected.
    assert int(result.loc[result["term"] == "study:subject", "df"].iloc[0]) > 0


def test_adonis2_rejects_a_zero_df_term_in_sequential_mode() -> None:
    """Sequential fitting has no such excuse -- the term is simply redundant."""
    beta, metadata = confounded_fixture()
    with pytest.raises(MicrobiomeSuiteError, match="aliased"):
        adonis2(beta, metadata, formula="d ~ study:subject + study", permutations=0)


def test_adonis2_marks_an_aliased_term_rather_than_leaving_a_blank_row() -> None:
    """A zero-df row is a real answer, but it must not read like an absent one."""
    beta, metadata = confounded_fixture()
    metadata = metadata.copy()
    metadata["study_copy"] = metadata["study"]
    result = adonis2(beta, metadata, formula="d ~ study + study_copy", permutations=0, by="margin")
    for term in ("study", "study_copy"):
        row = result.loc[result["term"] == term].iloc[0]
        assert int(row["df"]) == 0
        assert "aliased" in row["note"]
    assert result.loc[result["term"] == "Residual", "note"].iloc[0] == ""

    # The model itself still fits: the residual matches the non-degenerate formula.
    reference = adonis2(beta, metadata, formula="d ~ study", permutations=0)
    assert float(result.loc[result["term"] == "Residual", "r2"].iloc[0]) == pytest.approx(
        float(reference.loc[reference["term"] == "Residual", "r2"].iloc[0])
    )


@pytest.mark.parametrize("value", [None, 3, ["margin"]])
def test_adonis2_rejects_a_non_string_by(value: object) -> None:
    beta, metadata = fixture_beta_and_metadata()
    with pytest.raises(MicrobiomeSuiteError, match="terms.*margin"):
        adonis2(beta, metadata, formula="d ~ body_site", permutations=0, by=value)


def _repeated_measures_fixture(sizes: list[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Subjects of the given sizes, split into two treatment groups by subject.

    Treatment is constant within a subject (the between-subject case restricted
    permutation exists for) and carries a strong, deterministic separation, so
    the observed pseudo-F is large and any failure to restrict permutation
    properly is easy to see in the resulting p-value.
    """
    rows = []
    for subject_index, size in enumerate(sizes):
        treatment = "A" if subject_index % 2 == 0 else "B"
        offset = 5.0 if treatment == "A" else -5.0
        for replicate in range(size):
            rows.append(
                {
                    "sample_id": f"s{subject_index}-{replicate}",
                    "subject": f"subject-{subject_index}",
                    "treatment": treatment,
                    "feature": offset + 0.01 * replicate,
                }
            )
    metadata = pd.DataFrame(rows).set_index("sample_id")
    values = metadata["feature"].to_numpy()
    distances = np.abs(values[:, None] - values[None, :])
    beta = pd.DataFrame(distances, index=metadata.index, columns=metadata.index)
    return beta, metadata


def test_adonis2_rejects_a_between_subject_design_with_no_achievable_exchange() -> None:
    """Every plot a distinct size means no between-plot swap is possible; p=1.0 would lie."""
    beta, metadata = _repeated_measures_fixture([1, 2, 3, 4, 5, 6])
    with pytest.raises(MicrobiomeSuiteError, match="distinct size"):
        adonis2(
            beta,
            metadata,
            formula="d ~ treatment + (1 | subject)",
            permutations=99,
            seed=1,
        )


def test_adonis2_between_subject_design_with_repeated_sizes_still_works() -> None:
    beta, metadata = _repeated_measures_fixture([2, 2, 3, 3, 4, 4])
    result = adonis2(
        beta,
        metadata,
        formula="d ~ treatment + (1 | subject)",
        permutations=99,
        seed=1,
    )
    treatment_row = result.loc[result["term"] == "treatment"].iloc[0]
    assert float(treatment_row["p_value"]) < 1.0
    assert int(treatment_row["achievable_plot_permutations"]) == 8


def test_adonis2_warns_when_achievable_permutations_are_coarser_than_requested() -> None:
    beta, metadata = _repeated_measures_fixture([2, 2, 3, 3])
    with pytest.warns(UserWarning, match="4 distinct"):
        result = adonis2(
            beta,
            metadata,
            formula="d ~ treatment + (1 | subject)",
            permutations=99,
            seed=1,
        )
    treatment_row = result.loc[result["term"] == "treatment"].iloc[0]
    assert int(treatment_row["achievable_plot_permutations"]) == 4


def test_achievable_plot_permutations_matches_a_hand_checked_layout() -> None:
    """Sizes 2, 2, 3, 3 -> one size-2 pair and one size-3 pair -> 2! * 2! = 4."""
    scheme = PermutationScheme(plots=np.array([0, 0, 1, 1, 2, 2, 2, 3, 3, 3]))
    assert _achievable_plot_permutations(scheme, 10, cap=1000) == 4
