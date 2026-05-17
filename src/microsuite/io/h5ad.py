from __future__ import annotations

from pathlib import Path
from typing import Any, cast

import anndata as ad


def read_h5ad(path: Path) -> ad.AnnData:
    return ad.read_h5ad(path)


def write_h5ad(adata: ad.AnnData, path: Path) -> None:
    cast(Any, adata).write_h5ad(str(path))
