from __future__ import annotations

from typing import Any, cast

import anndata as ad
import numpy as np
from scipy import sparse


def dense_counts(adata: ad.AnnData) -> np.ndarray:
    matrix = adata.X
    if sparse.issparse(matrix):
        matrix = cast(Any, matrix).toarray()
    return np.asarray(matrix, dtype=np.float64)
