from __future__ import annotations

from pathlib import Path

from microsuite.io.metadata import read_indexed_tsv
from microsuite.io.taxonomy import normalize_taxonomy_columns


def test_qiime_metadata_type_row_is_ignored(tmp_path: Path) -> None:
    path = tmp_path / "metadata.tsv"
    path.write_text(
        "sample-id\tbody-site\n#q2:types\tcategorical\nsample-1\tgut\n",
        encoding="utf-8",
    )

    result = read_indexed_tsv(path, index_name="sample")

    assert result.index.tolist() == ["sample-1"]


def test_qiime_taxon_column_maps_to_taxonomy() -> None:
    mapping = normalize_taxonomy_columns(["Feature ID", "Taxon", "Confidence"])

    assert mapping == {"Taxon": "taxonomy"}
