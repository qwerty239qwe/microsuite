from __future__ import annotations

import shutil
import tempfile
from importlib.resources import files
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.diversity._matrix import dense_counts
from microsuite.io.h5ad import read_h5ad
from microsuite.runtime.runner import CommandLog, run_command

SUPPORTED_BACKENDS = ("native-correlation", "sparcc", "spieceasi", "flashweave")

SPIECEASI_SCRIPT = files("microsuite.networks.r").joinpath("spieceasi_network.R")
FLASHWEAVE_SCRIPT = files("microsuite.networks.julia").joinpath("flashweave_network.jl")


def network(
    *,
    backend: str,
    table: Path,
    output: Path,
    method: str = "pearson",
    transform: str = "relative",
    min_abs_weight: float = 0.3,
    min_prevalence: float = 0.1,
    top_n: int | None = None,
    pseudocount: float = 1.0,
    spieceasi_method: str = "mb",
    lambda_min_ratio: float = 0.01,
    nlambda: int = 20,
    sensitive: bool = False,
    heterogeneous: bool = False,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    if backend == "native-correlation":
        result = correlation_network(
            read_h5ad(ensure_input(table)),
            method=method,
            transform=transform,
            min_abs_weight=min_abs_weight,
            min_prevalence=min_prevalence,
            top_n=top_n,
            pseudocount=pseudocount,
        )
        result.to_csv(prepare_output(output, force=force), sep="\t", index=False)
        return
    if backend == "sparcc":
        result = sparcc_network(
            read_h5ad(ensure_input(table)),
            min_abs_weight=min_abs_weight,
            min_prevalence=min_prevalence,
            top_n=top_n,
            pseudocount=pseudocount,
        )
        result.to_csv(prepare_output(output, force=force), sep="\t", index=False)
        return
    if backend == "spieceasi":
        run_spieceasi(
            table=table,
            output=output,
            method=spieceasi_method,
            lambda_min_ratio=lambda_min_ratio,
            nlambda=nlambda,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    if backend == "flashweave":
        run_flashweave(
            table=table,
            output=output,
            sensitive=sensitive,
            heterogeneous=heterogeneous,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    raise MicrobiomeSuiteError(
        f"Unsupported network backend '{backend}'. Choose one of: {', '.join(SUPPORTED_BACKENDS)}"
    )


def correlation_network(
    adata: ad.AnnData,
    *,
    method: str = "pearson",
    transform: str = "relative",
    min_abs_weight: float = 0.3,
    min_prevalence: float = 0.1,
    top_n: int | None = None,
    pseudocount: float = 1.0,
) -> pd.DataFrame:
    matrix, feature_ids = _filtered_feature_matrix(
        adata,
        transform=transform,
        min_prevalence=min_prevalence,
        pseudocount=pseudocount,
    )
    return _edge_list_from_matrix(
        matrix,
        feature_ids,
        method=method,
        min_abs_weight=min_abs_weight,
        top_n=top_n,
        backend="native-correlation",
    )


def sparcc_network(
    adata: ad.AnnData,
    *,
    min_abs_weight: float = 0.3,
    min_prevalence: float = 0.1,
    top_n: int | None = None,
    pseudocount: float = 1.0,
) -> pd.DataFrame:
    matrix, feature_ids = _filtered_feature_matrix(
        adata,
        transform="clr",
        min_prevalence=min_prevalence,
        pseudocount=pseudocount,
    )
    edges = _edge_list_from_matrix(
        matrix,
        feature_ids,
        method="pearson",
        min_abs_weight=min_abs_weight,
        top_n=top_n,
        backend="sparcc",
    )
    edges["p_value"] = np.nan
    return edges


def run_spieceasi(
    *,
    table: Path,
    output: Path,
    method: str,
    lambda_min_ratio: float,
    nlambda: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if method not in {"mb", "glasso"}:
        raise MicrobiomeSuiteError("--spieceasi-method must be mb or glasso.")
    if not 0 < lambda_min_ratio <= 1:
        raise MicrobiomeSuiteError("--lambda-min-ratio must be greater than 0 and <= 1.")
    if nlambda < 1:
        raise MicrobiomeSuiteError("--nlambda must be at least 1.")
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise MicrobiomeSuiteError("SPIEC-EASI network inference requires Rscript and SpiecEasi.")
    prepare_output(output, force=force)
    with tempfile.TemporaryDirectory(prefix="microsuite-spieceasi-") as temp:
        counts = Path(temp) / "sample-feature-table.tsv"
        _write_sample_feature_table(read_h5ad(ensure_input(table)), counts)
        command = [
            rscript,
            str(SPIECEASI_SCRIPT),
            str(counts),
            str(output),
            method,
            str(lambda_min_ratio),
            str(nlambda),
        ]
        run_command(
            command,
            "SPIEC-EASI network inference failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="network",
                backend="spieceasi",
                inputs={"table": str(table)},
                outputs={"output": str(output)},
            ),
        )


def run_flashweave(
    *,
    table: Path,
    output: Path,
    sensitive: bool,
    heterogeneous: bool,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    julia = shutil.which("julia")
    if julia is None:
        raise MicrobiomeSuiteError("FlashWeave network inference requires Julia with FlashWeave.")
    prepare_output(output, force=force)
    with tempfile.TemporaryDirectory(prefix="microsuite-flashweave-") as temp:
        counts = Path(temp) / "sample-feature-table.tsv"
        _write_sample_feature_table(read_h5ad(ensure_input(table)), counts)
        command = [
            julia,
            str(FLASHWEAVE_SCRIPT),
            str(counts),
            str(output),
            str(sensitive).lower(),
            str(heterogeneous).lower(),
        ]
        run_command(
            command,
            "FlashWeave network inference failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="network",
                backend="flashweave",
                inputs={"table": str(table)},
                outputs={"output": str(output)},
            ),
        )


def _filtered_feature_matrix(
    adata: ad.AnnData,
    *,
    transform: str,
    min_prevalence: float,
    pseudocount: float,
) -> tuple[np.ndarray, list[str]]:
    if not 0 <= min_prevalence <= 1:
        raise MicrobiomeSuiteError("--min-prevalence must be between 0 and 1.")
    if pseudocount <= 0:
        raise MicrobiomeSuiteError("--pseudocount must be greater than 0.")
    counts = dense_counts(adata).astype(float)
    prevalence = (counts > 0).mean(axis=0)
    keep = prevalence >= min_prevalence
    if keep.sum() < 2:
        raise MicrobiomeSuiteError("At least two features pass the prevalence filter.")
    counts = counts[:, keep]
    feature_ids = adata.var_names.astype(str).to_numpy()[keep].tolist()
    return _transform_counts(counts, transform=transform, pseudocount=pseudocount), feature_ids


def _transform_counts(counts: np.ndarray, *, transform: str, pseudocount: float) -> np.ndarray:
    transform = transform.lower()
    if transform == "counts":
        return counts
    totals = counts.sum(axis=1, keepdims=True)
    relative = np.divide(counts, totals, out=np.zeros_like(counts), where=totals > 0)
    if transform == "relative":
        return relative
    if transform == "clr":
        logged = np.log(counts + pseudocount)
        return logged - logged.mean(axis=1, keepdims=True)
    raise MicrobiomeSuiteError("--transform must be counts, relative, or clr.")


def _edge_list_from_matrix(
    matrix: np.ndarray,
    feature_ids: list[str],
    *,
    method: str,
    min_abs_weight: float,
    top_n: int | None,
    backend: str,
) -> pd.DataFrame:
    method = method.lower()
    if method not in {"pearson", "spearman"}:
        raise MicrobiomeSuiteError("--method must be pearson or spearman.")
    if not 0 <= min_abs_weight <= 1:
        raise MicrobiomeSuiteError("--min-abs-weight must be between 0 and 1.")
    rows = []
    for i, source in enumerate(feature_ids):
        for j in range(i + 1, len(feature_ids)):
            target = feature_ids[j]
            if method == "pearson":
                result = stats.pearsonr(matrix[:, i], matrix[:, j])
            else:
                result = stats.spearmanr(matrix[:, i], matrix[:, j])
            weight = float(result.statistic)
            if np.isnan(weight) or abs(weight) < min_abs_weight:
                continue
            rows.append(
                {
                    "source": source,
                    "target": target,
                    "weight": weight,
                    "abs_weight": abs(weight),
                    "p_value": float(result.pvalue),
                    "method": method,
                    "backend": backend,
                }
            )
    edges = pd.DataFrame(rows)
    if edges.empty:
        return _empty_edges()
    edges = edges.sort_values(["abs_weight", "source", "target"], ascending=[False, True, True])
    if top_n is not None:
        if top_n < 1:
            raise MicrobiomeSuiteError("--top-n must be at least 1 when provided.")
        edges = edges.head(top_n)
    return edges.reset_index(drop=True)


def _empty_edges() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "source": pd.Series(dtype="string"),
            "target": pd.Series(dtype="string"),
            "weight": pd.Series(dtype="float64"),
            "abs_weight": pd.Series(dtype="float64"),
            "p_value": pd.Series(dtype="float64"),
            "method": pd.Series(dtype="string"),
            "backend": pd.Series(dtype="string"),
        }
    )


def _write_sample_feature_table(adata: ad.AnnData, output: Path) -> None:
    frame = pd.DataFrame(
        dense_counts(adata),
        index=adata.obs_names.astype(str),
        columns=adata.var_names.astype(str),
    )
    frame.index.name = "sample_id"
    frame.to_csv(output, sep="\t")
