from __future__ import annotations

from pathlib import Path

import anndata as ad

from microsuite._errors import MicrobiomeSuiteError
from microsuite.io.h5ad import read_h5ad, write_h5ad


def read_table(path: str | Path) -> ad.AnnData:
    table_path = Path(path)
    if table_path.suffix.lower() != ".h5ad":
        raise MicrobiomeSuiteError("Python SDK read_table currently expects an .h5ad file.")
    return read_h5ad(table_path)


def write_table(adata: ad.AnnData, path: str | Path) -> None:
    table_path = Path(path)
    if table_path.suffix.lower() != ".h5ad":
        raise MicrobiomeSuiteError("Python SDK write_table currently expects an .h5ad file.")
    write_h5ad(adata, table_path)
