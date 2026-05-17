from __future__ import annotations

from pathlib import Path

import pytest

from microsuite._errors import MicrobiomeSuiteError
from microsuite.methods.diversity_calc import diversity_calc
from microsuite.methods.tax_classify import tax_classify


def test_tax_classify_qiime2_requires_classifier(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.qza"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--classifier"):
        tax_classify(backend="qiime2", rep_seqs=rep_seqs, output=tmp_path / "taxonomy.qza")


def test_tax_classify_planned_method_message(tmp_path: Path) -> None:
    rep_seqs = tmp_path / "rep-seqs.fastq"
    rep_seqs.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="registered but not implemented"):
        tax_classify(backend="kraken2", rep_seqs=rep_seqs, output=tmp_path / "taxonomy.tsv")


def test_diversity_calc_qiime2_requires_phylogeny_for_unifrac(tmp_path: Path) -> None:
    table = tmp_path / "table.qza"
    table.write_text("placeholder", encoding="utf-8")

    with pytest.raises(MicrobiomeSuiteError, match="--phylogeny"):
        diversity_calc(
            backend="qiime2",
            metric="weighted-unifrac",
            table=table,
            output=tmp_path / "weighted-unifrac.qza",
        )


def test_diversity_calc_qiime2_reports_missing_qiime(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    table = tmp_path / "table.qza"
    table.write_text("placeholder", encoding="utf-8")

    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MicrobiomeSuiteError, match="qiime"):
        diversity_calc(
            backend="qiime2",
            metric="bray-curtis",
            table=table,
            output=tmp_path / "bray-curtis.qza",
        )
