"""Marshal an AnnData through an R batch-correction backend and back.

The read-back is the delicate half. Backends may drop features and return rows
in their own order, so the corrected matrix is realigned by feature label, never
by position: a positional rebuild yields a complete table with the right labels
on the wrong numbers.
"""

from __future__ import annotations

import json
from functools import partial
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, cast

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.backends import BatchBackend, resolve_backend
from microsuite.batch.value_type import record_batch_correction
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.container import resolve_batch_image
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog


def run_batch_correction(
    adata: ad.AnnData,
    *,
    backend: str,
    batch: str,
    covariates: list[str] | None = None,
    target: str | None = None,
    extra_params: dict[str, Any] | None = None,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
) -> ad.AnnData:
    record = resolve_backend(backend, covariates=covariates, target=target)
    covariate_list = list(covariates or [])
    _validate_design(adata, batch=batch, covariates=covariate_list, target=target)

    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        metadata_path = temp / "metadata.tsv"
        params_path = temp / "params.json"
        corrected_path = temp / "corrected.tsv"

        pd.DataFrame(dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")

        params: dict[str, Any] = {
            "batch": batch,
            "covariates": covariate_list,
            "target": target,
        }
        params.update(extra_params or {})
        params_path.write_text(json.dumps(params, indent=2, sort_keys=True), encoding="utf-8")

        invoke_r_script(
            backend=record.name,
            script_package="microsuite.batch.r",
            script_name=record.script,
            resolve_image=partial(resolve_batch_image, image=record.image),
            positional=[counts_path, metadata_path, params_path, corrected_path],
            runtime=runtime,
            image=image,
            engine=engine,
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="batch_correct",
                backend=record.name,
                inputs={"batch": batch, "covariates": ",".join(covariate_list)},
                outputs={"output": str(corrected_path)},
            ),
            local_missing_message=(
                f"batch correct --backend {record.name} requires external Rscript and the "
                f"R package '{record.package}'. Install R, then {record.install_hint}, or "
                f"use --runtime docker with the {record.image} image."
            ),
        )

        corrected = _read_corrected(corrected_path)

    return _rebuild(
        adata,
        corrected,
        record=record,
        batch=batch,
        covariates=covariate_list,
        target=target,
    )


def _validate_design(
    adata: ad.AnnData, *, batch: str, covariates: list[str], target: str | None
) -> None:
    available = list(adata.obs.columns)
    for label, column in [("--batch-col", batch), ("--target-col", target)]:
        if column is not None and column not in available:
            raise MicrobiomeSuiteError(
                f"{label} '{column}' not found in sample metadata. Available: "
                f"{', '.join(map(str, available))}"
            )
    missing = [name for name in covariates if name not in available]
    if missing:
        raise MicrobiomeSuiteError(
            f"--covariates not found in sample metadata: {', '.join(missing)}. "
            f"Available: {', '.join(map(str, available))}"
        )

    for label, column in [("--batch-col", batch), *[("--covariates", c) for c in covariates]]:
        _reject_na(adata, column=column, label=label)
    if target is not None:
        _reject_na(adata, column=target, label="--target-col")

    batch_values = adata.obs[batch].astype(str)
    if batch_values.nunique() < 2:
        raise MicrobiomeSuiteError(
            f"'{batch}' has one batch level ({batch_values.iloc[0]}); there is nothing to correct."
        )
    for name in covariates:
        values = adata.obs[name].astype(str)
        crosstab = pd.crosstab(batch_values, values)
        # Perfectly confounded: each batch sees exactly one covariate level.
        if (crosstab > 0).sum(axis=1).max() == 1:
            raise MicrobiomeSuiteError(
                f"Covariate '{name}' is perfectly confounded with batch '{batch}': every "
                f"batch contains a single '{name}' level, so no model can separate their "
                f"effects. Drop the covariate, or correct a different grouping."
            )


def _reject_na(adata: ad.AnnData, *, column: str, label: str) -> None:
    values = adata.obs[column]
    na_mask = values.isna()
    if na_mask.any():
        offenders = list(map(str, adata.obs_names[na_mask][:5]))
        raise MicrobiomeSuiteError(
            f"{label} '{column}' has missing (NA) values for samples: {', '.join(offenders)}. "
            f"Batch correction requires a fully populated column."
        )


def _read_corrected(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        raise MicrobiomeSuiteError("The batch-correction backend produced no output table.")
    frame = pd.read_csv(path, sep="\t", index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    nan_mask = cast("pd.Series", frame.isna().any(axis=1))
    if nan_mask.any():
        affected = list(frame.index[nan_mask][:5])
        raise MicrobiomeSuiteError(
            f"The batch-correction backend returned NA values for features: {', '.join(affected)}"
        )
    return frame


def _rebuild(
    adata: ad.AnnData,
    corrected: pd.DataFrame,
    *,
    record: BatchBackend,
    batch: str,
    covariates: list[str],
    target: str | None,
) -> ad.AnnData:
    unknown = [name for name in corrected.index if name not in set(map(str, adata.var_names))]
    if unknown:
        raise MicrobiomeSuiteError(
            f"The backend returned features absent from the input table: {', '.join(unknown[:5])}"
        )
    missing_samples = [
        name for name in map(str, adata.obs_names) if name not in set(corrected.columns)
    ]
    if missing_samples:
        raise MicrobiomeSuiteError(
            f"The backend dropped samples from the corrected table: "
            f"{', '.join(missing_samples[:5])}"
        )

    kept = [name for name in map(str, adata.var_names) if name in set(corrected.index)]
    aligned = corrected.loc[kept, [str(name) for name in adata.obs_names]]

    values = aligned.to_numpy(dtype=float)
    if record.value_type == "counts" and not np.allclose(values, np.round(values), atol=1e-6):
        raise MicrobiomeSuiteError(
            f"'{record.name}' declares value_type='counts' but returned non-integer values. "
            f"Refusing to stamp fractional data as counts."
        )

    result = cast(Any, adata[:, kept]).copy()
    result.X = values.T
    record_batch_correction(
        result,
        value_type=record.value_type,
        backend=record.name,
        batch=batch,
        covariates=covariates,
        target=target,
    )
    return result
