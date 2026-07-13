from __future__ import annotations

import json
import shutil
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.runner import CommandLog, run_command

ANCOMBC_SCRIPT = files("microsuite.diffab.r").joinpath("ancombc.R")


def run_ancombc(
    adata: ad.AnnData,
    *,
    output: Path,
    group: str | None = None,
    fix_formula: str | None = None,
    rand_formula: str | None = None,
    reference: dict[str, str] | None = None,
    prv_cut: float = 0.10,
    lib_cut: int = 0,
    struc_zero: bool = False,
    neg_lb: bool = False,
    p_adj_method: str = "BH",
    global_test: bool = False,
    pairwise: bool = False,
    trend: bool = False,
    dunnet: bool = False,
    pseudo_sens: bool = True,
    n_cl: int = 1,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    resolved_fix = fix_formula or group
    if not resolved_fix:
        raise MicrobiomeSuiteError("Provide --fix-formula or --group for ANCOM-BC.")
    reference = reference or {}
    obs_cols = set(adata.obs.columns)
    referenced = ([group] if group else []) + list(reference)
    missing = [c for c in referenced if c not in obs_cols]
    if missing:
        raise MicrobiomeSuiteError(f"Metadata columns not found in obs: {missing}")

    rscript = shutil.which("Rscript")
    if rscript is None:
        raise MicrobiomeSuiteError(
            "ANCOM-BC requires external Rscript and the R packages 'ANCOMBC' and "
            "'jsonlite'. Install R and those packages, then rerun this command."
        )

    with TemporaryDirectory() as temp_dir:
        temp = Path(temp_dir)
        counts_path = temp / "counts.tsv"
        metadata_path = temp / "metadata.tsv"
        params_path = temp / "params.json"

        pd.DataFrame(dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")

        params = {
            "fix_formula": resolved_fix,
            "rand_formula": rand_formula,
            "group": group,
            "reference": reference,
            "prv_cut": prv_cut,
            "lib_cut": lib_cut,
            "struc_zero": struc_zero,
            "neg_lb": neg_lb,
            "p_adj_method": p_adj_method,
            "global": global_test,
            "pairwise": pairwise,
            "trend": trend,
            "dunnet": dunnet,
            "pseudo_sens": pseudo_sens,
            "n_cl": n_cl,
        }
        params_path.write_text(json.dumps(params, indent=2), encoding="utf-8")

        run_command(
            [
                rscript,
                str(ANCOMBC_SCRIPT),
                str(counts_path),
                str(metadata_path),
                str(params_path),
                str(output),
            ],
            failure_message="ANCOM-BC failed.",
            run_dir=run_dir,
            log=CommandLog(
                task="diff_abundance",
                backend="ancombc",
                inputs={"fix_formula": resolved_fix, "rand_formula": rand_formula or ""},
                outputs={"output": str(output)},
            ),
            timeout=timeout,
        )
