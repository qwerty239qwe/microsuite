from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile

from typer.testing import CliRunner

from microsuite.cli.app import app
from microsuite.io.h5ad import write_h5ad
from microsuite.io.tsv import read_tsv

FIXTURE = Path(__file__).parent / "fixtures" / "moving_pictures_small"


def test_cli_help() -> None:
    result = CliRunner().invoke(app, ["--help"])

    assert result.exit_code == 0
    assert "Unified microbiome" in result.stdout


def test_cli_tiny_pipeline(tmp_path: Path) -> None:
    runner = CliRunner()
    h5ad = tmp_path / "table.h5ad"
    alpha = tmp_path / "alpha.tsv"
    beta = tmp_path / "beta.tsv"
    pcoa = tmp_path / "pcoa.tsv"
    barplot = tmp_path / "barplot.png"

    result = runner.invoke(
        app,
        [
            "import",
            "tsv",
            str(FIXTURE / "table.tsv"),
            "--metadata",
            str(FIXTURE / "metadata.tsv"),
            "--taxonomy",
            str(FIXTURE / "taxonomy.tsv"),
            "-o",
            str(h5ad),
        ],
    )
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(
        app, ["diversity", "alpha", str(h5ad), "--metric", "shannon", "-o", str(alpha)]
    )
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(
        app, ["diversity", "beta", str(h5ad), "--metric", "bray-curtis", "-o", str(beta)]
    )
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(app, ["ordination", "pcoa", str(beta), "-o", str(pcoa)])
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(
        app, ["viz", "barplot", str(h5ad), "--level", "genus", "-o", str(barplot)]
    )
    assert result.exit_code == 0, result.stdout

    assert alpha.exists()
    assert beta.exists()
    assert pcoa.exists()
    assert barplot.exists()


def test_cli_qiime_inspect(tmp_path: Path) -> None:
    artifact = tmp_path / "table.qza"
    with ZipFile(artifact, "w") as archive:
        archive.writestr("uuid/VERSION", "QIIME 2\narchive: 5\nframework: 2024.10.0\n")
        archive.writestr(
            "uuid/metadata.yaml",
            "uuid: uuid\ntype: FeatureTable[Frequency]\nformat: BIOMV210DirFmt\n",
        )
        archive.write(FIXTURE / "table.tsv", "uuid/data/feature-table.tsv")

    result = CliRunner().invoke(app, ["qiime", "inspect", str(artifact)])

    assert result.exit_code == 0, result.stdout
    assert "FeatureTable[Frequency]" in result.stdout
    assert "feature-table.tsv" in result.stdout


def test_cli_workflow_list() -> None:
    result = CliRunner().invoke(app, ["workflow", "list"])

    assert result.exit_code == 0, result.stdout
    assert "table-summary" in result.stdout
    assert "moving-pictures" in result.stdout


def test_cli_workflow_moving_pictures(tmp_path: Path) -> None:
    output = tmp_path / "run"

    result = CliRunner().invoke(
        app,
        ["workflow", "moving-pictures", "--out", str(output)],
    )

    assert result.exit_code == 0, result.stdout
    assert (output / "table.h5ad").exists()
    assert (output / "run.json").exists()


def test_cli_method_catalog() -> None:
    result = CliRunner().invoke(app, ["methods"])

    assert result.exit_code == 0, result.stdout
    assert "tax_classify" in result.stdout
    assert "diversity_calc" in result.stdout
    assert "qiime2" in result.stdout


def test_cli_tax_classify_kraken2_requires_database(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.fastq"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    result = CliRunner().invoke(
        app,
        [
            "tax_classify",
            "--backend",
            "kraken2",
            "--rep-seqs",
            str(rep_seqs),
            "-o",
            str(tmp_path / "taxonomy.tsv"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
    assert "--classifier" in str(result.exception)


def test_cli_diversity_calc_requires_phylogeny_for_unifrac(tmp_path: Path) -> None:
    table = tmp_path / "table.qza"
    table.write_text("placeholder", encoding="utf-8")

    result = CliRunner().invoke(
        app,
        [
            "diversity_calc",
            "--backend",
            "qiime2",
            "--metric",
            "weighted-unifrac",
            "--table",
            str(table),
            "-o",
            str(tmp_path / "weighted-unifrac.qza"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
    assert "--phylogeny" in str(result.exception)


def test_cli_method_alias_still_works(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.fastq"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    result = CliRunner().invoke(
        app,
        [
            "tax_classify",
            "--method",
            "kraken2",
            "--rep-seqs",
            str(rep_seqs),
            "-o",
            str(tmp_path / "taxonomy.tsv"),
        ],
    )

    assert result.exit_code == 1
    assert result.exception is not None
    assert "--classifier" in str(result.exception)


def test_cli_table_ecology_commands(tmp_path: Path) -> None:
    adata = read_tsv(FIXTURE / "table.tsv", FIXTURE / "metadata.tsv", FIXTURE / "taxonomy.tsv")
    table = tmp_path / "table.h5ad"
    write_h5ad(adata, table)

    relative = tmp_path / "relative.h5ad"
    abundance = tmp_path / "abundance.tsv"
    shared = tmp_path / "shared.tsv"
    rarefied = tmp_path / "rarefied.h5ad"
    runner = CliRunner()

    result = runner.invoke(
        app,
        [
            "normalize",
            "--backend",
            "native",
            "--method",
            "relative",
            "--table",
            str(table),
            "-o",
            str(relative),
        ],
    )
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(
        app,
        [
            "abundance",
            "--backend",
            "native",
            "--table",
            str(table),
            "--level",
            "genus",
            "-o",
            str(abundance),
        ],
    )
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(
        app,
        [
            "shared_taxa",
            "--backend",
            "native",
            "--table",
            str(table),
            "--level",
            "genus",
            "--group",
            "body_site",
            "-o",
            str(shared),
        ],
    )
    assert result.exit_code == 0, result.stdout

    result = runner.invoke(
        app,
        [
            "rarefy",
            "--backend",
            "native",
            "--table",
            str(table),
            "--depth",
            "10",
            "-o",
            str(rarefied),
        ],
    )
    assert result.exit_code == 0, result.stdout

    assert relative.exists()
    assert abundance.exists()
    assert shared.exists()
    assert rarefied.exists()
