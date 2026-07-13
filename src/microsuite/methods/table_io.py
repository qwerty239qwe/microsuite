from __future__ import annotations

from pathlib import Path

import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity._matrix import dense_counts
from microsuite.io.h5ad import read_h5ad
from microsuite.io.tsv import read_matrix_tsv
from microsuite.methods.normalize import normalize_native


def _matrix_frame(matrix, var_names, obs_names) -> pd.DataFrame:
    frame = pd.DataFrame(
        matrix.T,
        index=pd.Index([str(v) for v in var_names], name="feature_id"),
        columns=pd.Index([str(s) for s in obs_names]),
    )
    return frame


def export_table(
    *,
    table: Path,
    output: Path,
    layer: str | None = None,
    metadata: Path | None = None,
    force: bool = False,
) -> None:
    adata = read_h5ad(ensure_input(table))
    if layer is not None:
        if layer not in adata.layers:
            available = ", ".join(adata.layers.keys()) or "(none)"
            raise MicrobiomeSuiteError(f"Layer '{layer}' not found; available layers: {available}.")
        matrix = adata.layers[layer]
    else:
        matrix = dense_counts(adata)
    frame = _matrix_frame(matrix, adata.var_names, adata.obs_names)
    frame.to_csv(prepare_output(output, force=force), sep="\t")
    if metadata is not None:
        obs = pd.DataFrame(adata.obs)
        obs.index = obs.index.astype(str)
        obs.index.name = "sample"
        obs.to_csv(prepare_output(metadata, force=force), sep="\t")


def normalize_table(
    *,
    method: str,
    input_path: Path,
    output: Path,
    target_sum: float = 1_000_000.0,
    pseudocount: float = 1.0,
    min_prevalence: float = 0.1,
    force: bool = False,
) -> None:
    adata = read_matrix_tsv(ensure_input(input_path))
    result = normalize_native(
        adata,
        method=method,
        target_sum=target_sum,
        pseudocount=pseudocount,
        min_prevalence=min_prevalence,
    )
    frame = _matrix_frame(dense_counts(result), result.var_names, result.obs_names)
    frame.to_csv(prepare_output(output, force=force), sep="\t")
