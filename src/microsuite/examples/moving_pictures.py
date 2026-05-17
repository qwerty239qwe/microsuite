from __future__ import annotations

from pathlib import Path

from microsuite.data.moving_pictures import copy_small_fixture
from microsuite.workflows.table_summary import run_table_summary


def run_example(output: Path, *, force: bool = False) -> None:
    output.mkdir(parents=True, exist_ok=True)
    data_dir = output / "data"
    copy_small_fixture(data_dir, force=force)

    table = data_dir / "table.tsv"
    metadata = data_dir / "metadata.tsv"
    taxonomy = data_dir / "taxonomy.tsv"
    run_table_summary(
        output=output,
        table=table,
        metadata=metadata,
        taxonomy=taxonomy,
        input_format="tsv",
        force=force,
    )
