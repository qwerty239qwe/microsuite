"""R/vegan implementations of formula-based beta-diversity tests.

The native beta-significance implementation intentionally remains in
``ecology.py``.  This module owns the external-runtime boundary for vegan so
that formulas and restricted permutations are delegated to vegan rather than
reimplemented in Python.
"""

from __future__ import annotations

import json
import re
import shutil
import warnings
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from microsuite._errors import MicrobiomeSuiteError
from microsuite.runtime.container import (
    PathMapper,
    build_container_command,
    host_user_spec,
    require_engine,
    resolve_ecology_image,
    resolve_image_digest,
)
from microsuite.runtime.runner import CommandLog, run_command

VEGAN_METHODS = ("adonis2", "anosim2")
_FORMULA_IDENTIFIER = re.compile(r"(?<![`A-Za-z0-9_.])([A-Za-z][A-Za-z0-9_.]*)")
_R_FORMULA_WORDS = {
    "FALSE",
    "I",
    "Inf",
    "NA",
    "NaN",
    "NULL",
    "TRUE",
    "log",
    "log10",
    "sqrt",
}
_RANDOM_EFFECT_TERM = re.compile(r"\([^()]*\|[^()]*\)")


def vegan_beta_significance(
    distance_matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    *,
    column: str | None = None,
    method: str = "adonis2",
    formula: str | None = None,
    strata: str | None = None,
    permutations: int = 999,
    seed: int = 0,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
    run_dir: Path | None = None,
    timeout: float | None = None,
    sidecar_dir: Path | None = None,
) -> pd.DataFrame:
    """Run vegan's ``adonis2`` or the ``anosim2`` compatibility entry point.

    vegan calls its ANOSIM implementation ``anosim`` (there is no upstream
    ``anosim2`` function).  microsuite exposes ``anosim2`` as the stable name
    for this backend and dispatches it to ``vegan::anosim``.

    For ``adonis2``, ``formula`` is the right-hand side of the R formula and
    may be combined with ``strata`` and ``permutations``.  ``strata`` accepts
    one metadata column or a colon-separated metadata interaction, such as
    ``batch:subject``.  It restricts permutation exchangeability; it is not a
    mixed-effects random intercept.
    """

    method = _normalise_method(method)
    if runtime not in {"local", "docker"}:
        raise MicrobiomeSuiteError(
            f"Unsupported --runtime '{runtime}' for vegan; choose 'local' or 'docker'."
        )
    if permutations < 0:
        raise MicrobiomeSuiteError("--permutations must be non-negative.")

    matrix, meta, dropped_metadata_only = _align_inputs(distance_matrix, metadata)
    formula_text = _resolve_formula(column, formula)
    if method == "adonis2" and _RANDOM_EFFECT_TERM.search(formula_text):
        raise MicrobiomeSuiteError(
            "vegan adonis2 does not support lme4-style random-effects terms such as "
            "'(1 | batch:subject)'. Use strata='batch:subject' to restrict permutations; "
            "strata does not fit a random-effects model."
        )
    formula_columns = _formula_columns(formula_text, meta.columns)
    missing = sorted(name for name in formula_columns if name not in meta.columns)
    if missing:
        raise MicrobiomeSuiteError(
            f"Formula column(s) not found in sample metadata: {', '.join(missing)}"
        )
    strata_columns = _strata_columns(strata, meta.columns)

    used_columns = list(dict.fromkeys([*formula_columns, *strata_columns]))
    if used_columns:
        values = meta.loc[:, used_columns]
        if values.isna().any().any():
            missing_rows = values.index[values.isna().any(axis=1)].astype(str).tolist()
            raise MicrobiomeSuiteError(
                "Formula/strata metadata contains missing values for sample(s): "
                + ", ".join(missing_rows)
            )

    meta, strata_name = _materialise_strata(meta, strata, strata_columns)
    effective_permutations = permutations
    if strata_name is not None and permutations > 0:
        effective_permutations = _validate_strata_permutations(
            meta[strata_name], requested=permutations
        )
        _warn_formula_columns_constant_within_strata(
            meta,
            strata_name=strata_name,
            formula_columns=formula_columns,
        )

    if method == "anosim2":
        group_column = _anosim_group_column(formula_text, column, meta.columns)
        groups = meta[group_column].astype(str)
        if groups.nunique() < 2:
            raise MicrobiomeSuiteError("ANOSIM requires at least two groups.")
        if groups.value_counts().min() < 2:
            raise MicrobiomeSuiteError("ANOSIM requires at least two samples in every group.")

    with TemporaryDirectory(prefix="microsuite-vegan-") as temp_dir:
        temp = Path(temp_dir)
        distance_path = temp / "distance.tsv"
        metadata_path = temp / "metadata.tsv"
        result_path = temp / "result.tsv"
        matrix.to_csv(distance_path, sep="\t", index=True, index_label="sample_id")
        meta.to_csv(metadata_path, sep="\t", index=True, index_label="sample_id")
        _invoke_vegan(
            method=method,
            distance_path=distance_path,
            metadata_path=metadata_path,
            formula=formula_text,
            strata=strata_name,
            permutations=permutations,
            seed=seed,
            output_path=result_path,
            runtime=runtime,
            image=image,
            engine=engine,
            run_dir=run_dir,
            timeout=timeout,
            sidecar_dir=sidecar_dir,
        )
        result = pd.read_csv(result_path, sep="\t")

    if result.empty:
        raise MicrobiomeSuiteError(f"vegan {method} returned no result rows.")
    result["strata"] = strata if strata is not None else np.nan
    result["requested_permutations"] = permutations
    result["effective_permutations"] = effective_permutations
    result["n_samples"] = len(meta)
    result["dropped_metadata_only"] = ",".join(dropped_metadata_only)
    return result


def _normalise_method(method: str) -> str:
    value = method.lower().replace("-", "")
    aliases = {"adonis": "adonis2", "permanova": "adonis2", "anosim": "anosim2"}
    value = aliases.get(value, value)
    if value not in VEGAN_METHODS:
        raise MicrobiomeSuiteError(
            "vegan --method must be adonis2 or anosim2 (anosim is accepted as an alias)."
        )
    return value


def _resolve_formula(column: str | None, formula: str | None) -> str:
    if formula is None or not formula.strip():
        if column is None or not column.strip():
            raise MicrobiomeSuiteError("vegan requires --formula or --column.")
        formula = column
    formula = formula.strip()
    if "\n" in formula or "\r" in formula or "~" in formula:
        raise MicrobiomeSuiteError("--formula must be the right-hand side of an R formula.")
    return formula


def _strata_columns(strata: str | None, metadata_columns: pd.Index) -> list[str]:
    if strata is None:
        return []
    text = strata.strip()
    known = {str(column) for column in metadata_columns}
    if text in known:
        return [text]
    parts = [part.strip().strip("`") for part in text.split(":")]
    if len(parts) > 1 and all(part in known for part in parts):
        return list(dict.fromkeys(parts))
    raise MicrobiomeSuiteError(
        "Strata must be a metadata column or a colon-separated interaction of "
        f"metadata columns: {strata}"
    )


def _materialise_strata(
    metadata: pd.DataFrame, strata: str | None, strata_columns: list[str]
) -> tuple[pd.DataFrame, str | None]:
    if strata is None or len(strata_columns) <= 1:
        return metadata, strata_columns[0] if strata_columns else None
    result = metadata.copy()
    derived = "__microsuite_strata__"
    while derived in result.columns:
        derived = f"_{derived}"
    result[derived] = result.loc[:, strata_columns].astype(str).agg("\x1f".join, axis=1)
    return result, derived


def _validate_strata_permutations(strata: pd.Series, *, requested: int) -> int:
    """Validate block-restricted exchange and return the usable permutation count.

    vegan treats ``strata`` as non-movable blocks and freely permutes samples
    within each block.  For block sizes n_i, the complete permutation space is
    ``prod(factorial(n_i))`` including the observed ordering.  The calculation
    is capped at ``requested`` because only the number vegan can actually use
    is part of the public result.
    """
    counts = strata.value_counts(sort=False, dropna=False)
    singleton_count = int((counts == 1).sum())
    movable_count = int((counts > 1).sum())
    if movable_count == 0:
        raise MicrobiomeSuiteError(
            "Restricted permutation is impossible: every strata level contains one sample, "
            "so there are no legal non-identity permutations. Add repeated observations per "
            "stratum, choose a coarser --strata definition, or use unrestricted permutation."
        )

    if singleton_count:
        warnings.warn(
            f"{singleton_count} of {len(counts)} strata level(s) contain one sample and cannot "
            "move under restricted permutation.",
            UserWarning,
            stacklevel=2,
        )

    # Stop multiplying once at least requested + identity arrangements exist.
    arrangement_cap = requested + 1
    arrangements = 1
    for size in counts.astype(int):
        for factor in range(2, size + 1):
            arrangements *= factor
            if arrangements >= arrangement_cap:
                return requested

    effective = arrangements - 1
    if effective < requested:
        warnings.warn(
            "Restricted permutation can generate only "
            f"{effective} legal non-identity permutation(s), fewer than the requested "
            f"{requested}; vegan will use the complete set and p-value resolution is limited "
            f"to 1/{effective + 1} steps.",
            UserWarning,
            stacklevel=2,
        )
    return effective


def _warn_formula_columns_constant_within_strata(
    metadata: pd.DataFrame,
    *,
    strata_name: str,
    formula_columns: list[str],
) -> None:
    constant: list[str] = []
    grouped = metadata.groupby(strata_name, sort=False, dropna=False)
    for column in formula_columns:
        if bool((grouped[column].nunique(dropna=False) <= 1).all()):
            constant.append(column)
    if constant:
        warnings.warn(
            "Formula metadata column(s) are constant within every permutation stratum: "
            f"{', '.join(constant)}. Terms depending only on these columns have no "
            "within-stratum exchangeability, so their permutation p-values may be "
            "uninformative; sums of squares and R-squared remain descriptive.",
            UserWarning,
            stacklevel=2,
        )


def _formula_columns(formula: str, metadata_columns: pd.Index) -> list[str]:
    names = re.findall(r"`([^`]+)`", formula)
    for match in _FORMULA_IDENTIFIER.finditer(formula):
        token = match.group(1)
        end = match.end(1)
        if token in _R_FORMULA_WORDS or (formula[end:].lstrip().startswith("(")):
            continue
        names.append(token)
    known = set(str(column) for column in metadata_columns)
    return list(
        dict.fromkeys(name for name in names if name in known or name not in _R_FORMULA_WORDS)
    )


def _anosim_group_column(formula: str, column: str | None, metadata_columns: pd.Index) -> str:
    if column is not None:
        if formula != column and formula.strip("`") != column:
            raise MicrobiomeSuiteError(
                "ANOSIM supports one grouping column, not a multi-term formula."
            )
        return column
    names = [
        name for name in _formula_columns(formula, metadata_columns) if name in metadata_columns
    ]
    if len(names) != 1 or formula.strip("`") != names[0]:
        raise MicrobiomeSuiteError(
            "ANOSIM requires one grouping column: pass --column or a formula containing "
            "exactly one metadata column."
        )
    return names[0]


def _align_inputs(
    distance_matrix: pd.DataFrame, metadata: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    matrix = distance_matrix.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    if matrix.index.has_duplicates or matrix.columns.has_duplicates:
        raise MicrobiomeSuiteError("Distance matrix sample IDs must be unique.")
    if matrix.shape[0] != matrix.shape[1] or set(matrix.index) != set(matrix.columns):
        raise MicrobiomeSuiteError("Distance matrix must be square with matching sample IDs.")
    matrix = matrix.loc[matrix.index, matrix.index]
    values = matrix.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise MicrobiomeSuiteError("Distance matrix must contain only finite values.")
    if not np.allclose(values, values.T):
        raise MicrobiomeSuiteError("Distance matrix must be symmetric.")

    meta = metadata.copy()
    meta.index = meta.index.astype(str)
    if meta.index.has_duplicates:
        raise MicrobiomeSuiteError("Sample metadata IDs must be unique.")
    known_metadata = set(meta.index)
    missing_samples = [sample for sample in matrix.index if sample not in known_metadata]
    if missing_samples:
        shown = ", ".join(missing_samples[:10])
        suffix = "" if len(missing_samples) <= 10 else f", and {len(missing_samples) - 10} more"
        raise MicrobiomeSuiteError(
            f"Distance matrix contains {len(missing_samples)} sample ID(s) missing from "
            f"metadata: {shown}{suffix}"
        )
    sample_ids = list(matrix.index)
    if len(sample_ids) < 3:
        raise MicrobiomeSuiteError(
            "vegan beta-significance requires at least three shared samples."
        )
    dropped = [sample for sample in meta.index if sample not in set(matrix.index)]
    return matrix.loc[sample_ids, sample_ids], meta.loc[sample_ids], dropped


def _invoke_vegan(
    *,
    method: str,
    distance_path: Path,
    metadata_path: Path,
    formula: str,
    strata: str | None,
    permutations: int,
    seed: int,
    output_path: Path,
    runtime: str,
    image: str | None,
    engine: str,
    run_dir: Path | None,
    timeout: float | None,
    sidecar_dir: Path | None,
) -> None:
    positional: list[str | Path] = [
        distance_path,
        metadata_path,
        method,
        formula,
        strata or "",
        str(permutations),
        str(seed),
        output_path,
    ]
    log = CommandLog(
        task="beta_significance",
        backend=f"vegan-{method}",
        inputs={"formula": formula, "strata": strata},
        outputs={"output": str(output_path)},
        params={"permutations": permutations, "seed": seed, "runtime": runtime},
    )
    local_missing_message = (
        "vegan beta-significance requires external Rscript and the R package 'vegan'. "
        "Install R and vegan, then rerun this command, or use --runtime docker with "
        "the r-ecology image."
    )

    if runtime == "local":
        rscript = shutil.which("Rscript")
        if rscript is None:
            raise MicrobiomeSuiteError(local_missing_message)
        script = files("microsuite.diversity.r").joinpath("beta_significance.R")
        command = [rscript, str(script), *[str(arg) for arg in positional]]
        run_command(
            command,
            failure_message=f"vegan {method} failed.",
            run_dir=run_dir,
            log=log,
            timeout=timeout,
        )
        return

    resolved_image = resolve_ecology_image(image)
    require_engine(engine)
    paths = [arg for arg in positional if isinstance(arg, Path)]
    mapper = PathMapper()
    mountpoints: dict[Path, str] = {}

    def _mount(host_dir: Path, mode: str) -> None:
        resolved = host_dir.resolve()
        if resolved not in mountpoints:
            mountpoints[resolved] = f"/mnt/d{len(mountpoints)}"
        mapper.add_dir(host_dir, mode, mountpoints[resolved])

    for path in paths[:-1]:
        _mount(path.parent, "ro")
    _mount(paths[-1].parent, "rw")
    inner = ["/opt/microsuite/beta_significance.R"]
    inner.extend(mapper.to_container(arg) if isinstance(arg, Path) else arg for arg in positional)
    command = build_container_command(
        inner,
        resolved_image,
        mapper.mounts(),
        engine=engine,
        user=host_user_spec(),
    )
    run_command(
        command,
        failure_message=f"vegan {method} failed.",
        run_dir=run_dir,
        log=log,
        timeout=timeout,
    )
    if sidecar_dir is not None:
        sidecar_dir.mkdir(parents=True, exist_ok=True)
        (sidecar_dir / "vegan_container.json").write_text(
            json.dumps(
                {
                    "runtime": "docker",
                    "engine": engine,
                    "image": resolved_image,
                    "digest": resolve_image_digest(engine, resolved_image),
                    "method": method,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
