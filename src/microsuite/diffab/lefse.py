from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.batch.value_type import require_value_types
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.container import resolve_diffab_image
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog

P_ADJUST_METHODS = ("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
RESULT_COLUMNS = ("features", "scores")


def _validate_probability(value: float, name: str) -> None:
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise MicrobiomeSuiteError(f"LEfSe {name} must be between 0 and 1.")


def _validate_design(
    adata: ad.AnnData,
    *,
    group: str,
    subclass: str | None,
    reference: str | None,
) -> tuple[pd.DataFrame, str, str]:
    if group not in adata.obs.columns:
        raise MicrobiomeSuiteError(f"Group column not found in sample metadata: {group}")
    if subclass is not None and subclass not in adata.obs.columns:
        raise MicrobiomeSuiteError(f"Subclass column not found in sample metadata: {subclass}")
    if not adata.obs_names.is_unique:
        raise MicrobiomeSuiteError("LEfSe requires unique sample names.")
    if not adata.var_names.is_unique:
        raise MicrobiomeSuiteError("LEfSe requires unique feature names.")

    metadata = pd.DataFrame(adata.obs).copy()
    if metadata[group].isna().any():
        raise MicrobiomeSuiteError(f"LEfSe group column '{group}' contains missing values.")
    metadata[group] = metadata[group].astype(str)
    if metadata[group].str.strip().eq("").any():
        raise MicrobiomeSuiteError(f"LEfSe group column '{group}' contains empty values.")
    levels = sorted(metadata[group].unique().tolist())
    if len(levels) != 2:
        raise MicrobiomeSuiteError(
            f"LEfSe requires exactly two groups; '{group}' contains {len(levels)}."
        )
    counts = metadata[group].value_counts()
    if (counts < 2).any():
        raise MicrobiomeSuiteError("LEfSe requires at least two samples in each group.")

    if reference is None:
        resolved_reference = levels[0]
    else:
        resolved_reference = str(reference)
        if resolved_reference not in levels:
            raise MicrobiomeSuiteError(
                f"LEfSe reference '{resolved_reference}' is not a level of '{group}': "
                f"{', '.join(levels)}."
            )
    comparison = next(level for level in levels if level != resolved_reference)

    if subclass is not None:
        if metadata[subclass].isna().any():
            raise MicrobiomeSuiteError(
                f"LEfSe subclass column '{subclass}' contains missing values."
            )
        metadata[subclass] = metadata[subclass].astype(str)
        if metadata[subclass].str.strip().eq("").any():
            raise MicrobiomeSuiteError(f"LEfSe subclass column '{subclass}' contains empty values.")
        if metadata[subclass].nunique() < 2:
            raise MicrobiomeSuiteError("LEfSe subclass must contain at least two levels.")
        combinations = pd.crosstab(metadata[subclass], metadata[group])
        if combinations.shape[1] != 2 or (combinations == 0).any(axis=None):
            raise MicrobiomeSuiteError(
                "LEfSe subclass levels must be represented in both groups; subclass is a "
                "crossed blocking/replicate factor, not a nested or random-effect term."
            )

    return metadata, resolved_reference, comparison


def _validate_matrix(adata: ad.AnnData) -> np.ndarray:
    require_value_types(adata, ("counts", "relative"), operation="diff_abundance --backend lefse")
    matrix = dense_counts(adata)
    if matrix.ndim != 2 or matrix.shape != adata.shape:
        raise MicrobiomeSuiteError("LEfSe input must be a two-dimensional feature table.")
    if not np.isfinite(matrix).all():
        raise MicrobiomeSuiteError("LEfSe input contains non-finite abundance values.")
    if (matrix < 0).any():
        raise MicrobiomeSuiteError(
            "LEfSe requires non-negative counts or relative abundances; negative values "
            "such as CLR coordinates are not valid input."
        )
    zero_samples = np.asarray(matrix.sum(axis=1) <= 0).ravel()
    if zero_samples.any():
        names = ", ".join(map(str, adata.obs_names[zero_samples][:5]))
        raise MicrobiomeSuiteError(f"LEfSe samples have zero total abundance: {names}.")
    return matrix


def _params_sidecar(output: Path) -> Path:
    return output.with_name(f"{output.name}.params.json")


def run_lefse(
    adata: ad.AnnData,
    *,
    output: Path,
    group: str,
    subclass: str | None = None,
    reference: str | None = None,
    seed: int = 1234,
    kruskal_threshold: float = 0.05,
    wilcoxon_threshold: float = 0.05,
    lda_threshold: float = 2.0,
    p_adjust_method: str = "none",
    trim_names: bool = False,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
) -> None:
    """Run lefser with deterministic class ordering and validated inputs."""
    metadata, resolved_reference, comparison = _validate_design(
        adata, group=group, subclass=subclass, reference=reference
    )
    matrix = _validate_matrix(adata)
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed <= 2_147_483_647:
        raise MicrobiomeSuiteError("LEfSe seed must be an integer between 0 and 2147483647.")
    _validate_probability(kruskal_threshold, "kruskal_threshold")
    _validate_probability(wilcoxon_threshold, "wilcoxon_threshold")
    if not math.isfinite(lda_threshold) or lda_threshold < 0:
        raise MicrobiomeSuiteError("LEfSe lda_threshold must be non-negative.")
    method_lookup: dict[str, str] = {method.lower(): method for method in P_ADJUST_METHODS}
    try:
        resolved_method = method_lookup[p_adjust_method.lower()]
    except (AttributeError, KeyError):
        raise MicrobiomeSuiteError(
            f"Unsupported LEfSe p_adjust_method '{p_adjust_method}'; choose one of: "
            f"{', '.join(P_ADJUST_METHODS)}."
        ) from None

    output = Path(output)
    params_sidecar = _params_sidecar(output)
    existing = [path for path in (output, params_sidecar) if path.exists()]
    if existing and not force:
        raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {existing[0]}")
    output.parent.mkdir(parents=True, exist_ok=True)

    params = {
        "group": group,
        "subclass": subclass,
        "reference": resolved_reference,
        "comparison": comparison,
        "seed": seed,
        "kruskal_threshold": kruskal_threshold,
        "wilcoxon_threshold": wilcoxon_threshold,
        "lda_threshold": lda_threshold,
        "p_adjust_method": resolved_method,
        "trim_names": trim_names,
    }

    with (
        TemporaryDirectory() as input_temp_dir,
        TemporaryDirectory(dir=output.parent, prefix=".microsuite-lefse-") as stage_temp_dir,
    ):
        input_temp = Path(input_temp_dir)
        stage_temp = Path(stage_temp_dir)
        counts_path = input_temp / "counts.tsv"
        metadata_path = input_temp / "metadata.tsv"
        params_path = input_temp / "params.json"
        staged_output = stage_temp / "result.tsv"

        pd.DataFrame(matrix.T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )
        metadata.to_csv(metadata_path, sep="\t")
        params_path.write_text(
            json.dumps(params, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )

        invoke_r_script(
            backend="lefse",
            script_package="microsuite.diffab.r",
            script_name="lefse",
            resolve_image=resolve_diffab_image,
            positional=[counts_path, metadata_path, params_path, staged_output],
            runtime=runtime,
            image=image,
            engine=engine,
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="diff_abundance",
                backend="lefse",
                inputs=params,
                outputs={"output": str(output), "params": str(params_sidecar)},
            ),
            local_missing_message=(
                "LEfSe requires external Rscript and the R packages 'lefser' and "
                "'jsonlite'. Install R and those packages, then rerun this command, or "
                "use --runtime docker with the r-diffab-lefse image."
            ),
        )

        if not staged_output.is_file() or staged_output.stat().st_size == 0:
            raise MicrobiomeSuiteError("LEfSe did not produce its required result table.")
        try:
            result = pd.read_csv(staged_output, sep="\t")
        except Exception as exc:
            raise MicrobiomeSuiteError(f"LEfSe produced an unreadable result table: {exc}") from exc
        if tuple(result.columns) != RESULT_COLUMNS:
            raise MicrobiomeSuiteError(
                "LEfSe result schema changed: expected columns "
                f"{', '.join(RESULT_COLUMNS)}, found {', '.join(map(str, result.columns))}."
            )
        numeric_scores = pd.to_numeric(result["scores"], errors="coerce")
        if not result.empty and not np.isfinite(numeric_scores).all():
            raise MicrobiomeSuiteError("LEfSe produced non-numeric or non-finite LDA scores.")

        staged_params = stage_temp / params_sidecar.name
        staged_params.write_text(
            json.dumps(params, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        if force:
            for path in (output, params_sidecar):
                if path.exists():
                    if path.is_dir():
                        shutil.rmtree(path)
                    else:
                        path.unlink()
        staged_output.replace(output)
        staged_params.replace(params_sidecar)
        container_sidecar = stage_temp / "lefse_container.json"
        if container_sidecar.exists():
            container_sidecar.replace(output.parent / container_sidecar.name)
