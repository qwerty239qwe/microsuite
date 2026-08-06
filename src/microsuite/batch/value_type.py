"""The output-scale contract for batch-corrected tables.

Batch correction backends disagree about what they return: ComBat-seq and
ConQuR emit integer counts, MMUPHin and MetaDICT emit relative abundances, and
PLSDA-batch emits CLR log-ratios. Downstream, ANCOM-BC and ALDEx2 require
counts, ``rarefy`` requires counts, and ``normalize`` must not transform data
that is already transformed.

Nothing about a float matrix reveals which of the three it is, so the producing
command records it and the consuming commands assert on it. An absent key means
"unknown" and is always allowed: every table written before 0.3.0 lacks it, and
this contract must never change the behaviour of existing pipelines.
"""

from __future__ import annotations

import anndata as ad

from microsuite._errors import MicrobiomeSuiteError

VALUE_TYPES: tuple[str, ...] = ("counts", "relative", "clr")


def record_batch_correction(
    adata: ad.AnnData,
    *,
    value_type: str,
    backend: str,
    batch: str,
    covariates: list[str],
    target: str | None,
) -> None:
    """Stamp the corrected table with its scale and how it was produced."""
    if value_type not in VALUE_TYPES:
        raise MicrobiomeSuiteError(
            f"Unknown value_type '{value_type}'. Choose one of: {', '.join(VALUE_TYPES)}"
        )
    info = adata.uns.get("microsuite")
    if not isinstance(info, dict):
        info = {}
    info = dict(info)
    info["value_type"] = value_type
    info["batch_correct"] = {
        "backend": backend,
        "batch": batch,
        "covariates": list(covariates),
        "target": target,
    }
    adata.uns["microsuite"] = info


def read_value_type(adata: ad.AnnData) -> str | None:
    """Return the recorded scale, or None when the table does not declare one."""
    info = adata.uns.get("microsuite")
    if not isinstance(info, dict):
        return None
    value_type = info.get("value_type")
    return value_type if isinstance(value_type, str) else None


def require_value_types(adata: ad.AnnData, allowed: tuple[str, ...], *, operation: str) -> None:
    """Raise when the table declares a scale that ``operation`` cannot consume."""
    value_type = read_value_type(adata)
    if value_type is None or value_type in allowed:
        return
    info = adata.uns.get("microsuite")
    provenance = info.get("batch_correct") if isinstance(info, dict) else None
    origin = ""
    if isinstance(provenance, dict) and provenance.get("backend"):
        origin = (
            f" It was produced by 'batch correct --backend {provenance['backend']}', "
            f"which emits '{value_type}'."
        )
    raise MicrobiomeSuiteError(
        f"{operation} requires a table of type {' or '.join(allowed)}, "
        f"but this table is '{value_type}'.{origin} "
        f"Use a backend that accepts '{value_type}', or correct a different way."
    )
