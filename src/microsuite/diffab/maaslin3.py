from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

import anndata as ad
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.diversity._matrix import dense_counts
from microsuite.runtime.container import resolve_diffab_image
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog

ABUNDANCE_RESULTS = "abundance_results.tsv"
PREVALENCE_RESULTS = "prevalence_results.tsv"
NORMALIZATIONS = ("TSS", "CLR", "NONE")
TRANSFORMS = ("LOG", "PLOG", "NONE")


def _resolved_formula(
    *,
    group: str | None,
    formula: str | None,
    fix_formula: str | None,
    rand_formula: str | None,
) -> str:
    if formula and any(value for value in (group, fix_formula, rand_formula)):
        raise MicrobiomeSuiteError(
            "MaAsLin 3 --formula cannot be combined with --group, --fix-formula, "
            "or --rand-formula. Put the complete fixed/random design in --formula."
        )
    if group and fix_formula:
        raise MicrobiomeSuiteError("MaAsLin 3 accepts either --group or --fix-formula, not both.")

    if formula:
        resolved = formula.strip()
    else:
        fixed = (fix_formula or group or "").strip()
        if not fixed:
            raise MicrobiomeSuiteError(
                "Provide --formula, --fix-formula, or --group for MaAsLin 3."
            )
        random = (rand_formula or "").strip()
        resolved = f"{fixed} + {random}" if random else fixed

    if not resolved:
        raise MicrobiomeSuiteError("MaAsLin 3 formula cannot be empty.")
    return resolved if resolved.startswith("~") else f"~ {resolved}"


def _validate_options(
    *, normalization: str, transform: str, min_prevalence: float, min_abundance: float
) -> tuple[str, str]:
    normalization = normalization.upper()
    transform = transform.upper()
    if normalization not in NORMALIZATIONS:
        raise MicrobiomeSuiteError(
            f"Unsupported MaAsLin 3 normalization '{normalization}'; "
            f"choose one of: {', '.join(NORMALIZATIONS)}."
        )
    if transform not in TRANSFORMS:
        raise MicrobiomeSuiteError(
            f"Unsupported MaAsLin 3 transform '{transform}'; "
            f"choose one of: {', '.join(TRANSFORMS)}."
        )
    if not math.isfinite(min_prevalence) or not 0 <= min_prevalence <= 1:
        raise MicrobiomeSuiteError("MaAsLin 3 min_prevalence must be between 0 and 1.")
    if not math.isfinite(min_abundance) or min_abundance < 0:
        raise MicrobiomeSuiteError("MaAsLin 3 min_abundance must be non-negative.")
    return normalization, transform


def _normalize_reference(reference: str | None, *, metadata_columns: set[str]) -> str | None:
    if reference is None:
        return None
    if not reference.strip():
        raise MicrobiomeSuiteError(
            "MaAsLin 3 reference cannot be empty; use 'column,level' pairs separated by ';'."
        )

    normalized: list[str] = []
    seen: set[str] = set()
    for raw_pair in reference.split(";"):
        if raw_pair.count(",") != 1:
            raise MicrobiomeSuiteError(
                "Invalid MaAsLin 3 reference; expected 'column,level' pairs separated by ';'."
            )
        column, level = (part.strip() for part in raw_pair.split(",", maxsplit=1))
        if not column or not level:
            raise MicrobiomeSuiteError(
                "Invalid MaAsLin 3 reference; column and level must both be non-empty."
            )
        if column in seen:
            raise MicrobiomeSuiteError(
                f"MaAsLin 3 reference column is specified more than once: {column}"
            )
        if column not in metadata_columns:
            raise MicrobiomeSuiteError(
                f"MaAsLin 3 reference column not found in sample metadata: {column}"
            )
        seen.add(column)
        normalized.append(f"{column},{level}")
    return ";".join(normalized)


def _replace_output(staged: Path, output: Path, *, force: bool) -> None:
    if output.exists():
        if not force:
            raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {output}")
        if output.is_dir():
            shutil.rmtree(output)
        else:
            output.unlink()
    staged.replace(output)


def run_maaslin3(
    adata: ad.AnnData,
    *,
    output: Path,
    group: str | None = None,
    formula: str | None = None,
    fix_formula: str | None = None,
    rand_formula: str | None = None,
    normalization: str = "TSS",
    transform: str = "LOG",
    min_prevalence: float = 0.0,
    min_abundance: float = 0.0,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    reference: str | None = None,
    engine: str = "docker",
) -> None:
    """Run MaAsLin 3 and write its two model families to an output directory."""
    resolved_formula = _resolved_formula(
        group=group,
        formula=formula,
        fix_formula=fix_formula,
        rand_formula=rand_formula,
    )
    normalization, transform = _validate_options(
        normalization=normalization,
        transform=transform,
        min_prevalence=min_prevalence,
        min_abundance=min_abundance,
    )
    reference = _normalize_reference(
        reference,
        metadata_columns={str(column) for column in adata.obs.columns},
    )
    if group and group not in adata.obs.columns:
        raise MicrobiomeSuiteError(f"Group column not found in sample metadata: {group}")
    if output.exists() and not force:
        raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    params = {
        "formula": resolved_formula,
        **({"reference": reference} if reference is not None else {}),
        "normalization": normalization,
        "transform": transform,
        "min_prevalence": min_prevalence,
        "min_abundance": min_abundance,
    }

    with (
        TemporaryDirectory() as input_temp_dir,
        TemporaryDirectory(dir=output.parent, prefix=".microsuite-maaslin3-") as stage_temp_dir,
    ):
        input_temp = Path(input_temp_dir)
        stage_temp = Path(stage_temp_dir)
        counts_path = input_temp / "counts.tsv"
        metadata_path = input_temp / "metadata.tsv"
        params_path = input_temp / "params.json"
        staged_output = stage_temp / "result"

        pd.DataFrame(dense_counts(adata).T, index=adata.var_names, columns=adata.obs_names).to_csv(
            counts_path, sep="\t"
        )
        pd.DataFrame(adata.obs).to_csv(metadata_path, sep="\t")
        params_path.write_text(json.dumps(params, indent=2), encoding="utf-8")

        invoke_r_script(
            backend="maaslin3",
            script_package="microsuite.diffab.r",
            script_name="maaslin3",
            resolve_image=resolve_diffab_image,
            positional=[counts_path, metadata_path, params_path, staged_output],
            runtime=runtime,
            image=image,
            engine=engine,
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="diff_abundance",
                backend="maaslin3",
                inputs=params,
                outputs={
                    "abundance": str(output / ABUNDANCE_RESULTS),
                    "prevalence": str(output / PREVALENCE_RESULTS),
                },
            ),
            local_missing_message=(
                "MaAsLin 3 requires external Rscript and the R packages 'maaslin3' and "
                "'jsonlite'. Install R and those packages, then rerun this command, or use "
                "--runtime docker with the r-diffab-maaslin3 image."
            ),
        )

        for filename in (ABUNDANCE_RESULTS, PREVALENCE_RESULTS):
            result = staged_output / filename
            if not result.is_file() or result.stat().st_size == 0:
                raise MicrobiomeSuiteError(
                    f"MaAsLin 3 did not produce its required result table: {filename}"
                )

        sidecar = stage_temp / "maaslin3_container.json"
        if sidecar.exists():
            sidecar.replace(staged_output / sidecar.name)
        _replace_output(staged_output, output, force=force)
