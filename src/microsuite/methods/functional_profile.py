from __future__ import annotations

import csv
import gzip
import hashlib
import importlib
import json
import math
import os
import re
import shutil
import subprocess
import uuid
from collections.abc import Sequence
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.container import (
    PathMapper,
    build_container_command,
    host_user_spec,
    require_engine,
    resolve_functional_image,
    resolve_image_digest,
    resolve_picrust2_image,
)
from microsuite.runtime.r_backend import invoke_r_script
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("picrust2", "tax4fun2", "humann")

TAX4FUN2_SCRIPT = files("microsuite.functional.r").joinpath("tax4fun2.R")
TAX4FUN2_VERSION = "1.1.5"
TAX4FUN2_FUNCTIONS = "functional_prediction.tsv"
TAX4FUN2_PATHWAYS = "pathway_prediction.tsv"
TAX4FUN2_COVERAGE = "coverage.tsv"
TAX4FUN2_MANIFEST = "tax4fun2_manifest.json"
TAX4FUN2_REQUIRED_OUTPUTS = (
    TAX4FUN2_FUNCTIONS,
    TAX4FUN2_PATHWAYS,
    TAX4FUN2_COVERAGE,
    TAX4FUN2_MANIFEST,
)
PICRUSt2_DATABASES = ("SC", "oldIMG", "custom")
PICRUSt2_MANIFEST = "picrust2_manifest.json"
PICRUSt2_STANDARD_OUTPUTS = {
    "ec": "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    "ko": "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
    "weighted_nsti": "EC_metagenome_out/weighted_nsti.tsv.gz",
    "pathways": "pathways_out/path_abun_unstrat.tsv.gz",
    "coverage": "pathways_out/path_cov_unstrat.tsv.gz",
}


def functional_profile(
    *,
    backend: str,
    output_dir: Path,
    table: Path | None = None,
    rep_seqs: Path | None = None,
    reads: Path | None = None,
    database: Path | None = None,
    protein_database: Path | None = None,
    threads: int | str = "1",
    database_mode: str = "Ref99NR",
    min_identity: float = 0.97,
    normalize_pathways: bool = False,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
    picrust2_database: str = "SC",
    picrust2_ref_dir1: Path | None = None,
    picrust2_ref_dir2: Path | None = None,
    picrust2_custom_trait_tables: Path | str | Sequence[Path | str] | None = None,
    picrust2_custom_trait_tables_ref1: Path | str | Sequence[Path | str] | None = None,
    picrust2_custom_trait_tables_ref2: Path | str | Sequence[Path | str] | None = None,
    picrust2_marker_gene_table: Path | None = None,
    picrust2_marker_gene_table_ref1: Path | None = None,
    picrust2_marker_gene_table_ref2: Path | None = None,
    picrust2_pathway_map: Path | None = None,
    picrust2_reaction_func: Path | str | None = None,
    picrust2_regroup_map: Path | None = None,
    picrust2_no_regroup: bool = False,
    picrust2_no_pathways: bool = False,
    picrust2_coverage: bool = False,
    picrust2_max_nsti: float = 2.0,
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "functional_profile")
    if backend not in {"tax4fun2", "picrust2"} and (
        runtime != "local" or image is not None or engine != "docker"
    ):
        raise MicrobiomeSuiteError(
            "--runtime, --image, and --engine apply only to --backend tax4fun2 or picrust2."
        )
    if backend == "picrust2":
        functional_profile_picrust2(
            table=table,
            rep_seqs=rep_seqs,
            output_dir=output_dir,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
            runtime=runtime,
            image=image,
            engine=engine,
            picrust2_database=picrust2_database,
            picrust2_ref_dir1=picrust2_ref_dir1,
            picrust2_ref_dir2=picrust2_ref_dir2,
            picrust2_custom_trait_tables=picrust2_custom_trait_tables,
            picrust2_custom_trait_tables_ref1=picrust2_custom_trait_tables_ref1,
            picrust2_custom_trait_tables_ref2=picrust2_custom_trait_tables_ref2,
            picrust2_marker_gene_table=picrust2_marker_gene_table,
            picrust2_marker_gene_table_ref1=picrust2_marker_gene_table_ref1,
            picrust2_marker_gene_table_ref2=picrust2_marker_gene_table_ref2,
            picrust2_pathway_map=picrust2_pathway_map,
            picrust2_reaction_func=picrust2_reaction_func,
            picrust2_regroup_map=picrust2_regroup_map,
            picrust2_no_regroup=picrust2_no_regroup,
            picrust2_no_pathways=picrust2_no_pathways,
            picrust2_coverage=picrust2_coverage,
            picrust2_max_nsti=picrust2_max_nsti,
        )
        return
    if backend == "tax4fun2":
        functional_profile_tax4fun2(
            table=table,
            rep_seqs=rep_seqs,
            database=database,
            output_dir=output_dir,
            threads=threads,
            database_mode=database_mode,
            min_identity=min_identity,
            normalize_pathways=normalize_pathways,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
            runtime=runtime,
            image=image,
            engine=engine,
        )
        return
    if backend == "humann":
        functional_profile_humann(
            reads=reads,
            database=database,
            protein_database=protein_database,
            output_dir=output_dir,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return


def functional_profile_picrust2(
    *,
    table: Path | None,
    rep_seqs: Path | None,
    output_dir: Path,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
    runtime: str = "local",
    image: str | None = None,
    engine: str = "docker",
    picrust2_database: str = "SC",
    picrust2_ref_dir1: Path | None = None,
    picrust2_ref_dir2: Path | None = None,
    picrust2_custom_trait_tables: Path | str | Sequence[Path | str] | None = None,
    picrust2_custom_trait_tables_ref1: Path | str | Sequence[Path | str] | None = None,
    picrust2_custom_trait_tables_ref2: Path | str | Sequence[Path | str] | None = None,
    picrust2_marker_gene_table: Path | None = None,
    picrust2_marker_gene_table_ref1: Path | None = None,
    picrust2_marker_gene_table_ref2: Path | None = None,
    picrust2_pathway_map: Path | None = None,
    picrust2_reaction_func: Path | str | None = None,
    picrust2_regroup_map: Path | None = None,
    picrust2_no_regroup: bool = False,
    picrust2_no_pathways: bool = False,
    picrust2_coverage: bool = False,
    picrust2_max_nsti: float = 2.0,
) -> None:
    if table is None:
        raise MicrobiomeSuiteError("--table is required for --backend picrust2.")
    if rep_seqs is None:
        raise MicrobiomeSuiteError("--rep-seqs is required for --backend picrust2.")
    table = ensure_input(table).resolve()
    rep_seqs = ensure_input(rep_seqs).resolve()
    _validate_picrust2_inputs(table, rep_seqs)
    mode = _canonical_picrust2_database(picrust2_database)
    config = _validate_picrust2_config(
        mode=mode,
        ref_dir1=picrust2_ref_dir1,
        ref_dir2=picrust2_ref_dir2,
        trait_tables=picrust2_custom_trait_tables,
        trait_tables_ref1=picrust2_custom_trait_tables_ref1,
        trait_tables_ref2=picrust2_custom_trait_tables_ref2,
        marker_gene_table=picrust2_marker_gene_table,
        marker_gene_table_ref1=picrust2_marker_gene_table_ref1,
        marker_gene_table_ref2=picrust2_marker_gene_table_ref2,
        pathway_map=picrust2_pathway_map,
        reaction_func=picrust2_reaction_func,
        regroup_map=picrust2_regroup_map,
        no_regroup=picrust2_no_regroup,
    )
    if not math.isfinite(picrust2_max_nsti) or picrust2_max_nsti < 0:
        raise MicrobiomeSuiteError("--picrust2-max-nsti must be a finite non-negative number.")
    if picrust2_coverage and picrust2_no_pathways:
        raise MicrobiomeSuiteError(
            "--picrust2-coverage cannot be used with --picrust2-no-pathways because "
            "pathway coverage requires pathway inference."
        )
    resolved_threads = resolve_threads(threads)
    output_dir = Path(output_dir)
    _validate_output_target(output_dir, force=force)
    external_paths = [table, rep_seqs, *config["external_paths"]]
    _validate_picrust2_path_collisions(output_dir, external_paths)

    script_name = config["script_name"]
    executable_name = script_name
    if runtime not in {"local", "docker"}:
        raise MicrobiomeSuiteError(
            f"Unsupported --runtime '{runtime}' for picrust2; choose 'local' or 'docker'."
        )
    resolved_image: str | None = None
    image_digest: str | None = None
    executable: str | None = None
    if runtime == "local":
        executable = shutil.which(executable_name)
        if executable is None:
            raise MicrobiomeSuiteError(
                f"PICRUSt2 functional profiling requires '{executable_name}'. "
                "Install PICRUSt2 or use --runtime docker."
            )
    else:
        resolved_image = resolve_picrust2_image(image)
        if resolved_image is None:
            raise MicrobiomeSuiteError("PICRUSt2 Docker image resolution returned no image.")
        require_engine(engine)
    picrust2_version = _probe_picrust2_version(
        runtime=runtime,
        executable=executable,
        script_name=script_name,
        image=resolved_image,
        engine=engine,
        timeout=timeout,
    )
    if mode == "SC" and _picrust2_version_tuple(picrust2_version) < (2, 6):
        raise MicrobiomeSuiteError(
            "PICRUSt2-SC requires PICRUSt2 >= 2.6; the selected pipeline reports "
            f"{picrust2_version}."
        )

    params: dict[str, Any] = {
        "picrust2_database": mode,
        "picrust2_max_nsti": picrust2_max_nsti,
        "picrust2_coverage": picrust2_coverage,
        "picrust2_no_regroup": picrust2_no_regroup,
        "picrust2_no_pathways": picrust2_no_pathways,
        "threads": resolved_threads,
        "runtime": runtime,
    }
    if config["custom"]:
        params["picrust2_ref_dir1"] = str(config["ref_dir1"])
        if config["ref_dir2"] is not None:
            params["picrust2_ref_dir2"] = str(config["ref_dir2"])
    if config["pathway_map"] is not None:
        params["picrust2_pathway_map"] = str(config["pathway_map"])
    if config["reaction_func_path"] is not None:
        params["picrust2_reaction_func"] = str(config["reaction_func_path"])
    elif config["reaction_func_value"] is not None:
        params["picrust2_reaction_func"] = config["reaction_func_value"]
    if config["regroup_map"] is not None:
        params["picrust2_regroup_map"] = str(config["regroup_map"])
    outputs = _picrust2_log_outputs(
        output_dir, no_pathways=picrust2_no_pathways, coverage=picrust2_coverage
    )
    log_inputs: dict[str, Any] = {"table": str(table), "rep_seqs": str(rep_seqs)}
    if config["external_paths"]:
        log_inputs["database_inputs"] = [str(path) for path in config["external_paths"]]

    with TemporaryDirectory(
        dir=output_dir.parent, prefix=".microsuite-picrust2-"
    ) as stage_temp_dir:
        staged_output = Path(stage_temp_dir) / "result"
        staged_output.mkdir()
        if runtime == "local":
            assert executable is not None
            command = _picrust2_command(
                executable,
                table=table,
                rep_seqs=rep_seqs,
                output=staged_output,
                threads=resolved_threads,
                max_nsti=picrust2_max_nsti,
                coverage=picrust2_coverage,
                no_regroup=picrust2_no_regroup,
                no_pathways=picrust2_no_pathways,
                config=config,
            )
        else:
            mapper = _picrust2_mapper(
                table=table,
                rep_seqs=rep_seqs,
                output=staged_output,
                config=config,
            )
            inner = _picrust2_command(
                script_name,
                table=table,
                rep_seqs=rep_seqs,
                output=staged_output,
                threads=resolved_threads,
                max_nsti=picrust2_max_nsti,
                coverage=picrust2_coverage,
                no_regroup=picrust2_no_regroup,
                no_pathways=picrust2_no_pathways,
                config=config,
                mapper=mapper,
            )
            assert resolved_image is not None
            command = build_container_command(
                inner,
                resolved_image,
                mapper.mounts(),
                engine=engine,
                user=host_user_spec(),
                entrypoint="",
            )
            image_digest = resolve_image_digest(engine, resolved_image)
        run_command(
            command,
            "PICRUSt2 functional profiling failed.",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="functional_profile",
                backend="picrust2",
                inputs=log_inputs,
                outputs=outputs,
                params=params,
            ),
        )
        discovered = _validate_picrust2_outputs(
            staged_output,
            mode=mode,
            no_pathways=picrust2_no_pathways,
            coverage=picrust2_coverage,
        )
        manifest = _picrust2_manifest(
            table=table,
            rep_seqs=rep_seqs,
            mode=mode,
            params=params,
            discovered=discovered,
            config=config,
            runtime=runtime,
            engine=engine,
            image=resolved_image,
            digest=image_digest,
            picrust2_version=picrust2_version,
        )
        (staged_output / PICRUSt2_MANIFEST).write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        _replace_output_dir(staged_output, output_dir, force=force)


def _canonical_picrust2_database(value: str) -> str:
    if not isinstance(value, str):
        raise MicrobiomeSuiteError("--picrust2-database must be SC, oldIMG, or custom.")
    normalized = value.strip().lower()
    aliases = {"sc": "SC", "oldimg": "oldIMG", "custom": "custom"}
    if normalized not in aliases:
        raise MicrobiomeSuiteError(
            f"Unsupported --picrust2-database '{value}'. Choose SC, oldIMG, or custom."
        )
    return aliases[normalized]


_PICRUST2_VERSION_RE = re.compile(
    r"PICRUSt2(?:\s+version)?\s+v?([0-9]+\.[0-9]+(?:\.[0-9]+)?)", re.IGNORECASE
)
_GENERIC_VERSION_RE = re.compile(r"\bv?([0-9]+\.[0-9]+(?:\.[0-9]+)?)\b")


def _parse_picrust2_version(stdout: str, stderr: str) -> str:
    text = "\n".join(part for part in (stdout, stderr) if part)
    match = _PICRUST2_VERSION_RE.search(text) or _GENERIC_VERSION_RE.search(text)
    if match is None:
        raise MicrobiomeSuiteError(
            "PICRUSt2 version probe did not report a parseable version. "
            "Run the selected PICRUSt2 pipeline with --version to inspect it."
        )
    return match.group(1)


def _picrust2_version_tuple(version: str) -> tuple[int, int]:
    major, minor, *_ = (int(part) for part in version.split("."))
    return major, minor


def _probe_picrust2_version(
    *,
    runtime: str,
    executable: str | None,
    script_name: str,
    image: str | None,
    engine: str,
    timeout: float | None,
) -> str:
    if runtime == "local":
        if executable is None:
            raise MicrobiomeSuiteError("PICRUSt2 version probe has no local executable.")
        command = [executable, "--version"]
    else:
        if image is None:
            raise MicrobiomeSuiteError("PICRUSt2 version probe has no container image.")
        command = build_container_command(
            [script_name, "--version"],
            image,
            [],
            engine=engine,
            user=host_user_spec(),
            entrypoint="",
        )
    try:
        kwargs: dict[str, Any] = {
            "capture_output": True,
            "text": True,
            "check": False,
        }
        if timeout is not None:
            kwargs["timeout"] = timeout
        result = subprocess.run(command, **kwargs)
    except (OSError, subprocess.SubprocessError) as exc:
        raise MicrobiomeSuiteError(
            "PICRUSt2 version probe failed for the selected executable/image."
        ) from exc
    stdout = str(result.stdout or "")
    stderr = str(result.stderr or "")
    if result.returncode != 0:
        detail = stderr.strip() or stdout.strip() or "no diagnostic output"
        raise MicrobiomeSuiteError(
            f"PICRUSt2 version probe failed (exit {result.returncode}): {detail}"
        )
    return _parse_picrust2_version(stdout, stderr)


def _normalize_picrust2_paths(
    value: Path | str | Sequence[Path | str] | None, *, flag: str
) -> list[Path]:
    if value is None:
        return []
    values: list[Path | str]
    if isinstance(value, (Path, str)):
        values = [value]
    else:
        values = list(value)
    paths: list[Path] = []
    for item in values:
        for part in str(item).split(","):
            if part.strip():
                paths.append(Path(part.strip()).resolve())
    if len(paths) != len(set(paths)):
        raise MicrobiomeSuiteError(f"{flag} contains duplicate paths.")
    return paths


def _ensure_nonempty_file(path: Path, *, flag: str) -> Path:
    path = ensure_input(path).resolve()
    if path.stat().st_size == 0:
        raise MicrobiomeSuiteError(f"{flag} must be a non-empty file: {path}")
    return path


def _validate_picrust2_matrix_file(
    path: Path, *, flag: str
) -> tuple[Path, tuple[str, ...]]:
    """Validate an upstream custom trait/marker matrix without accepting junk files."""
    path = _ensure_nonempty_file(path, flag=flag)
    try:
        with _open_table_text(path) as handle:
            rows = csv.reader(handle, delimiter="\t")
            try:
                header = next(rows)
            except StopIteration as exc:
                raise MicrobiomeSuiteError(f"{flag} must contain a header and data rows.") from exc
            columns = [value.strip() for value in header]
            if len(columns) < 2 or any(not value for value in columns):
                raise MicrobiomeSuiteError(
                    f"{flag} must have a non-empty row-ID column and at least one trait column."
                )
            if len(columns) != len(set(columns)):
                raise MicrobiomeSuiteError(f"{flag} column IDs must be unique.")
            row_ids: set[str] = set()
            total = 0.0
            data_rows = 0
            for line_number, row in enumerate(rows, start=2):
                if not row or not any(value.strip() for value in row):
                    continue
                if len(row) != len(columns):
                    raise MicrobiomeSuiteError(
                        f"{flag} line {line_number} has {len(row)} fields; "
                        f"expected {len(columns)}."
                    )
                row_id = row[0].strip()
                if not row_id:
                    raise MicrobiomeSuiteError(f"{flag} has an empty row ID at line {line_number}.")
                if row_id in row_ids:
                    raise MicrobiomeSuiteError(f"{flag} row IDs must be unique: {row_id}")
                row_ids.add(row_id)
                data_rows += 1
                for value in row[1:]:
                    try:
                        number = float(value)
                    except ValueError as exc:
                        raise MicrobiomeSuiteError(
                            f"{flag} contains a non-numeric value at line {line_number}: {value!r}"
                        ) from exc
                    if not math.isfinite(number) or number < 0:
                        raise MicrobiomeSuiteError(
                            f"{flag} values must be finite and non-negative."
                        )
                    total += number
    except MicrobiomeSuiteError:
        raise
    except (OSError, UnicodeError, csv.Error) as exc:
        raise MicrobiomeSuiteError(f"{flag} is not a readable tab-delimited matrix.") from exc
    if not data_rows:
        raise MicrobiomeSuiteError(f"{flag} must contain at least one data row.")
    if total <= 0:
        raise MicrobiomeSuiteError(f"{flag} must have a nonzero numeric total.")
    return path, tuple(columns[1:])


def _picrust2_trait_family_id(path: Path) -> str:
    """Match PICRUSt2's one-splitext custom trait-family derivation exactly."""
    family_id = os.path.splitext(os.path.basename(str(path)))[0]
    if not family_id:
        raise MicrobiomeSuiteError(
            f"PICRUSt2 custom trait table has an empty family ID: {path.name}"
        )
    return family_id


def _validate_picrust2_reference(path: Path, *, flag: str) -> Path:
    path = _ensure_dir(path, flag=flag).resolve()
    basename = path.name
    fasta_candidates = [
        path / f"{basename}.fna.gz",
        path / f"{basename}.fasta.gz",
        path / f"{basename}.fna",
        path / f"{basename}.fasta",
    ]
    present = [candidate for candidate in fasta_candidates if candidate.is_file()]
    if len(present) != 1:
        raise MicrobiomeSuiteError(
            f"{flag} must contain exactly one of {basename}.fna.gz, {basename}.fasta.gz, "
            f"{basename}.fna, or {basename}.fasta; found {len(present)}."
        )
    required_files = [
        present[0],
        *(path / f"{basename}{suffix}" for suffix in (".tre", ".hmm", ".model")),
    ]
    missing = [file.name for file in required_files if not file.is_file()]
    empty = [file.name for file in required_files if file.is_file() and file.stat().st_size == 0]
    if missing:
        raise MicrobiomeSuiteError(
            f"{flag} is missing required EPA-ng reference files: " + ", ".join(missing)
        )
    if empty:
        raise MicrobiomeSuiteError(
            f"{flag} contains empty EPA-ng reference files: " + ", ".join(empty)
        )
    return path


def _validate_picrust2_config(
    *,
    mode: str,
    ref_dir1: Path | None,
    ref_dir2: Path | None,
    trait_tables: Path | str | Sequence[Path | str] | None,
    trait_tables_ref1: Path | str | Sequence[Path | str] | None,
    trait_tables_ref2: Path | str | Sequence[Path | str] | None,
    marker_gene_table: Path | None,
    marker_gene_table_ref1: Path | None,
    marker_gene_table_ref2: Path | None,
    pathway_map: Path | None,
    reaction_func: Path | str | None,
    regroup_map: Path | None,
    no_regroup: bool,
) -> dict[str, Any]:
    if no_regroup and regroup_map is not None:
        raise MicrobiomeSuiteError(
            "--picrust2-no-regroup cannot be combined with --picrust2-regroup-map."
        )
    if trait_tables is not None and trait_tables_ref1 is not None:
        raise MicrobiomeSuiteError(
            "Use either picrust2_custom_trait_tables or "
            "picrust2_custom_trait_tables_ref1, not both."
        )
    if marker_gene_table is not None and marker_gene_table_ref1 is not None:
        raise MicrobiomeSuiteError(
            "Use either picrust2_marker_gene_table or picrust2_marker_gene_table_ref1, not both."
        )
    traits1 = _normalize_picrust2_paths(
        trait_tables if trait_tables is not None else trait_tables_ref1,
        flag="--picrust2-custom-trait-tables-ref1",
    )
    traits2 = _normalize_picrust2_paths(
        trait_tables_ref2,
        flag="--picrust2-custom-trait-tables-ref2",
    )
    marker1 = marker_gene_table if marker_gene_table is not None else marker_gene_table_ref1
    custom_fields = (
        ref_dir1,
        ref_dir2,
        traits1,
        traits2,
        marker1,
        marker_gene_table_ref2,
    )
    has_custom_fields = any(value not in (None, [], ()) for value in custom_fields)
    if mode != "custom" and has_custom_fields:
        raise MicrobiomeSuiteError(
            f"Custom PICRUSt2 reference, trait, and marker options cannot be used with "
            f"--picrust2-database {mode}."
        )

    external_paths: list[Path] = []
    ref1 = ref2 = None
    if mode == "custom":
        if ref_dir1 is None:
            raise MicrobiomeSuiteError("--picrust2-ref-dir1 is required for custom PICRUSt2.")
        if not traits1:
            raise MicrobiomeSuiteError(
                "At least one --picrust2-custom-trait-tables-ref1 file is required "
                "for custom PICRUSt2."
            )
        if marker1 is None:
            raise MicrobiomeSuiteError(
                "--picrust2-marker-gene-table-ref1 is required for custom PICRUSt2."
            )
        has_ref2 = ref_dir2 is not None
        ref2_fields = (traits2, marker_gene_table_ref2)
        if has_ref2 != all(value not in (None, [], ()) for value in ref2_fields):
            raise MicrobiomeSuiteError(
                "Dual-reference custom PICRUSt2 requires ref_dir2, trait tables ref2, and "
                "marker_gene_table_ref2 together."
            )
        if not has_ref2 and any(value not in (None, [], ()) for value in ref2_fields):
            raise MicrobiomeSuiteError("ref2 custom PICRUSt2 options require picrust2_ref_dir2.")
        ref1 = _validate_picrust2_reference(ref_dir1, flag="--picrust2-ref-dir1")
        external_paths.append(ref1)
        if ref2 is not None or ref_dir2 is not None:
            ref2 = _validate_picrust2_reference(ref_dir2, flag="--picrust2-ref-dir2")
            external_paths.append(ref2)
        trait_results1 = [
            _validate_picrust2_matrix_file(
                path, flag="--picrust2-custom-trait-tables-ref1"
            )
            for path in traits1
        ]
        traits1 = [path for path, _ in trait_results1]
        trait_columns1: dict[str, tuple[str, ...]] = {}
        for path, columns in trait_results1:
            family_id = _picrust2_trait_family_id(path)
            if family_id in trait_columns1:
                raise MicrobiomeSuiteError(
                    "Duplicate PICRUSt2 custom trait family ID in reference 1: "
                    f"{family_id}"
                )
            trait_columns1[family_id] = columns
        marker1, marker_columns1 = _validate_picrust2_matrix_file(
            marker1, flag="--picrust2-marker-gene-table-ref1"
        )
        external_paths.extend(traits1)
        external_paths.append(marker1)
        if ref2 is not None:
            trait_results2 = [
                _validate_picrust2_matrix_file(
                    path, flag="--picrust2-custom-trait-tables-ref2"
                )
                for path in traits2
            ]
            traits2 = [path for path, _ in trait_results2]
            trait_columns2: dict[str, tuple[str, ...]] = {}
            for path, columns in trait_results2:
                family_id = _picrust2_trait_family_id(path)
                if family_id in trait_columns2:
                    raise MicrobiomeSuiteError(
                        "Duplicate PICRUSt2 custom trait family ID in reference 2: "
                        f"{family_id}"
                    )
                trait_columns2[family_id] = columns
            if marker_gene_table_ref2 is None:
                raise MicrobiomeSuiteError(
                    "--picrust2-marker-gene-table-ref2 is required for dual custom PICRUSt2."
                )
            marker_gene_table_ref2, marker_columns2 = _validate_picrust2_matrix_file(
                marker_gene_table_ref2, flag="--picrust2-marker-gene-table-ref2"
            )
            if set(trait_columns1) != set(trait_columns2):
                missing_ref2 = sorted(set(trait_columns1) - set(trait_columns2))
                missing_ref1 = sorted(set(trait_columns2) - set(trait_columns1))
                details: list[str] = []
                if missing_ref2:
                    details.append("missing from reference 2: " + ", ".join(missing_ref2))
                if missing_ref1:
                    details.append("missing from reference 1: " + ", ".join(missing_ref1))
                raise MicrobiomeSuiteError(
                    "Dual custom PICRUSt2 trait tables must have matching family IDs ("
                    + "; ".join(details)
                    + ")."
                )
            mismatched_family = next(
                (
                    family_id
                    for family_id in trait_columns1
                    if trait_columns1[family_id] != trait_columns2[family_id]
                ),
                None,
            )
            if mismatched_family is not None:
                raise MicrobiomeSuiteError(
                    "Dual custom PICRUSt2 trait tables must have matching column IDs for "
                    f"family {mismatched_family}."
                )
            if marker_columns1 != marker_columns2:
                raise MicrobiomeSuiteError(
                    "Dual custom PICRUSt2 marker tables must have matching column IDs."
                )
            external_paths.extend(traits2)
            external_paths.append(marker_gene_table_ref2)

    pathway_map_path = (
        _ensure_nonempty_file(pathway_map, flag="--picrust2-pathway-map")
        if pathway_map is not None
        else None
    )
    regroup_map_path = (
        _ensure_nonempty_file(regroup_map, flag="--picrust2-regroup-map")
        if regroup_map is not None
        else None
    )
    reaction_path: Path | None = None
    reaction_value: str | None = None
    if reaction_func is not None:
        if isinstance(reaction_func, Path):
            reaction_path = _ensure_nonempty_file(reaction_func, flag="--picrust2-reaction-func")
        else:
            reaction_text = str(reaction_func).strip()
            if not reaction_text:
                raise MicrobiomeSuiteError("--picrust2-reaction-func cannot be empty.")
            candidate = Path(reaction_text)
            looks_like_path = (
                candidate.exists()
                or "/" in reaction_text
                or "\\" in reaction_text
                or candidate.suffix in {".tsv", ".txt", ".gz", ".map"}
            )
            if looks_like_path:
                reaction_path = _ensure_nonempty_file(candidate, flag="--picrust2-reaction-func")
            else:
                reaction_value = reaction_text
    if pathway_map_path is not None:
        external_paths.append(pathway_map_path)
    if reaction_path is not None:
        external_paths.append(reaction_path)
    if regroup_map_path is not None:
        external_paths.append(regroup_map_path)

    dual = mode == "custom" and ref2 is not None
    return {
        "custom": mode == "custom",
        "dual": dual,
        "script_name": (
            "picrust2_pipeline.py" if mode == "SC" or dual else "picrust2_pipeline_singleRef.py"
        ),
        "ref_dir1": ref1,
        "ref_dir2": ref2,
        "trait_tables_ref1": traits1,
        "trait_tables_ref2": traits2,
        "marker_gene_table_ref1": marker1,
        "marker_gene_table_ref2": marker_gene_table_ref2,
        "pathway_map": pathway_map_path,
        "reaction_func_path": reaction_path,
        "reaction_func_value": reaction_value,
        "regroup_map": regroup_map_path,
        "external_paths": external_paths,
    }


def _picrust2_command(
    executable: str,
    *,
    table: Path,
    rep_seqs: Path,
    output: Path,
    threads: int,
    max_nsti: float,
    coverage: bool,
    no_regroup: bool,
    no_pathways: bool,
    config: dict[str, Any],
    mapper: PathMapper | None = None,
) -> list[str]:
    def path(value: Path) -> str:
        return mapper.to_container(value) if mapper is not None else str(value)

    output_arg = mapper.container_dir(output) if mapper is not None else str(output)
    command = [
        executable,
        "-s",
        path(rep_seqs),
        "-i",
        path(table),
        "-o",
        output_arg,
        "-p",
        str(threads),
    ]
    if config["custom"]:
        if config["dual"]:
            command.extend(
                [
                    "--ref_dir1",
                    mapper.container_dir(config["ref_dir1"]) if mapper else str(config["ref_dir1"]),
                ]
            )
            command.extend(
                [
                    "--ref_dir2",
                    mapper.container_dir(config["ref_dir2"]) if mapper else str(config["ref_dir2"]),
                ]
            )
            command.extend(
                [
                    "--custom_trait_tables_ref1",
                    ",".join(path(item) for item in config["trait_tables_ref1"]),
                ]
            )
            command.extend(
                [
                    "--custom_trait_tables_ref2",
                    ",".join(path(item) for item in config["trait_tables_ref2"]),
                ]
            )
            command.extend(["--marker_gene_table_ref1", path(config["marker_gene_table_ref1"])])
            command.extend(["--marker_gene_table_ref2", path(config["marker_gene_table_ref2"])])
        else:
            command.extend(
                [
                    "-r",
                    mapper.container_dir(config["ref_dir1"]) if mapper else str(config["ref_dir1"]),
                ]
            )
            command.extend(
                [
                    "--custom_trait_tables",
                    ",".join(path(item) for item in config["trait_tables_ref1"]),
                ]
            )
            command.extend(["--marker_gene_table", path(config["marker_gene_table_ref1"])])
    if config["pathway_map"] is not None:
        command.extend(["--pathway_map", path(config["pathway_map"])])
    if config["reaction_func_path"] is not None:
        command.extend(["--reaction_func", path(config["reaction_func_path"])])
    elif config["reaction_func_value"] is not None:
        command.extend(["--reaction_func", config["reaction_func_value"]])
    if config["regroup_map"] is not None:
        command.extend(["--regroup_map", path(config["regroup_map"])])
    if no_regroup:
        command.append("--no_regroup")
    if no_pathways:
        command.append("--no_pathways")
    if coverage:
        command.append("--coverage")
    command.extend(["--max_nsti", str(max_nsti)])
    return command


def _picrust2_mapper(
    *, table: Path, rep_seqs: Path, output: Path, config: dict[str, Any]
) -> PathMapper:
    mapper = PathMapper()
    if config["ref_dir1"] is not None:
        mapper.add_dir(config["ref_dir1"], "ro", "/microsuite/ref1")
    if config["ref_dir2"] is not None:
        mapper.add_dir(config["ref_dir2"], "ro", "/microsuite/ref2")
    files_to_mount = [
        table,
        rep_seqs,
        *config["trait_tables_ref1"],
        *config["trait_tables_ref2"],
        *[
            path
            for path in (config["marker_gene_table_ref1"], config["marker_gene_table_ref2"])
            if path
        ],
        *[
            path
            for path in (config["pathway_map"], config["reaction_func_path"], config["regroup_map"])
            if path
        ],
    ]
    for index, path in enumerate(dict.fromkeys(files_to_mount)):
        mapper.add_file(path, "ro", f"/microsuite/input{index}/{path.name}")
    mapper.add_dir(output, "rw", "/microsuite/output")
    return mapper


def _picrust2_log_outputs(output: Path, *, no_pathways: bool, coverage: bool) -> dict[str, str]:
    outputs = {"output_dir": str(output), "manifest": str(output / PICRUSt2_MANIFEST)}
    for name, relative in PICRUSt2_STANDARD_OUTPUTS.items():
        if name == "coverage" and not coverage:
            continue
        if name == "pathways" and no_pathways:
            continue
        outputs[name] = str(output / relative)
    return outputs


def _validate_picrust2_outputs(
    output: Path, *, mode: str, no_pathways: bool, coverage: bool
) -> list[str]:
    required = []
    if mode in {"SC", "oldIMG"}:
        required.extend(
            [
                PICRUSt2_STANDARD_OUTPUTS["ec"],
                PICRUSt2_STANDARD_OUTPUTS["ko"],
                PICRUSt2_STANDARD_OUTPUTS["weighted_nsti"],
            ]
        )
        if not no_pathways:
            required.append(PICRUSt2_STANDARD_OUTPUTS["pathways"])
    else:
        function_outputs = sorted(
            path.relative_to(output).as_posix()
            for path in output.glob("*_metagenome_out/pred_metagenome_unstrat.tsv.gz")
        )
        if not function_outputs:
            raise MicrobiomeSuiteError(
                "PICRUSt2 custom mode did not produce any *_metagenome_out/"
                "pred_metagenome_unstrat.tsv.gz function output."
            )
        required.extend(function_outputs)
    if coverage:
        required.append(PICRUSt2_STANDARD_OUTPUTS["coverage"])
    for relative in dict.fromkeys(required):
        _validate_gz_tsv(output / relative, relative)
    return sorted(
        path.relative_to(output).as_posix() for path in output.rglob("*") if path.is_file()
    )


def _validate_gz_tsv(path: Path, relative: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise MicrobiomeSuiteError(f"PICRUSt2 did not produce a required output: {relative}")
    try:
        with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
            header = handle.readline()
            row = handle.readline()
    except (OSError, EOFError, UnicodeError) as exc:
        raise MicrobiomeSuiteError(
            f"PICRUSt2 output is not a readable gzipped TSV: {relative}"
        ) from exc
    if "\t" not in header or not row.strip():
        raise MicrobiomeSuiteError(f"PICRUSt2 output is not a non-empty gzipped TSV: {relative}")


def _picrust2_manifest(
    *,
    table: Path,
    rep_seqs: Path,
    mode: str,
    params: dict[str, Any],
    discovered: list[str],
    config: dict[str, Any],
    runtime: str,
    engine: str,
    image: str | None,
    digest: str | None,
    picrust2_version: str,
) -> dict[str, Any]:
    custom_database: dict[str, Any] | None = None
    if config["custom"]:
        custom_database = {
            "ref_dir1": _fingerprint_path(config["ref_dir1"]),
            "ref_dir2": _fingerprint_path(config["ref_dir2"]) if config["ref_dir2"] else None,
            "trait_tables_ref1": [_fingerprint_path(path) for path in config["trait_tables_ref1"]],
            "trait_tables_ref2": [_fingerprint_path(path) for path in config["trait_tables_ref2"]],
            "marker_gene_table_ref1": _fingerprint_path(config["marker_gene_table_ref1"]),
            "marker_gene_table_ref2": (
                _fingerprint_path(config["marker_gene_table_ref2"])
                if config["marker_gene_table_ref2"]
                else None
            ),
            "pathway_map": _fingerprint_path(config["pathway_map"])
            if config["pathway_map"]
            else None,
            "reaction_func": _fingerprint_path(config["reaction_func_path"])
            if config["reaction_func_path"]
            else config["reaction_func_value"],
            "regroup_map": _fingerprint_path(config["regroup_map"])
            if config["regroup_map"]
            else None,
        }
    container = None
    if runtime == "docker":
        container = {"engine": engine, "image": image, "digest": digest}
    return {
        "schema_version": "microsuite-picrust2.v1",
        "picrust2_version": picrust2_version,
        "database_mode": mode,
        "parameters": params,
        "inputs": {"table": _fingerprint_path(table), "rep_seqs": _fingerprint_path(rep_seqs)},
        "custom_database": custom_database,
        "runtime": runtime,
        "container": container,
        "discovered_outputs": discovered,
    }


def _fingerprint_path(path: Path | None) -> dict[str, Any] | None:
    if path is None:
        return None
    path = Path(path).resolve()
    if path.is_file():
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return {
            "path": str(path),
            "kind": "file",
            "size": path.stat().st_size,
            "sha256": digest.hexdigest(),
        }
    files_data: dict[str, Any] = {}
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        fingerprint = _fingerprint_path(child)
        files_data[child.relative_to(path).as_posix()] = fingerprint
    return {"path": str(path), "kind": "directory", "files": files_data}


def _validate_output_target(path: Path, *, force: bool) -> None:
    if path.exists() and not path.is_dir():
        raise MicrobiomeSuiteError(f"Output path exists and is not a directory: {path}")
    if path.exists() and any(path.iterdir()) and not force:
        raise MicrobiomeSuiteError(f"Output directory exists, pass --force to overwrite: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def _validate_picrust2_path_collisions(output: Path, external_paths: Sequence[Path]) -> None:
    output = output.resolve()
    for external in external_paths:
        external = Path(external).resolve()
        if _path_contains(output, external) or _path_contains(external, output):
            raise MicrobiomeSuiteError(
                "PICRUSt2 output path collides with an input or database path: "
                f"{output} and {external}"
            )


def _path_contains(parent: Path, child: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def functional_profile_tax4fun2(
    *,
    table: Path | None,
    rep_seqs: Path | None,
    database: Path | None,
    output_dir: Path,
    threads: int | str,
    database_mode: str,
    min_identity: float,
    normalize_pathways: bool,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
    runtime: str,
    image: str | None,
    engine: str,
) -> None:
    if table is None:
        raise MicrobiomeSuiteError("--table is required for --backend tax4fun2.")
    if rep_seqs is None:
        raise MicrobiomeSuiteError("--rep-seqs is required for --backend tax4fun2.")
    if database is None:
        raise MicrobiomeSuiteError(
            "--database is required for --backend tax4fun2 and should point to "
            "Tax4Fun2_ReferenceData_v2."
        )
    if database_mode not in {"Ref99NR", "Ref100NR"}:
        raise MicrobiomeSuiteError("--database-mode must be Ref99NR or Ref100NR.")
    if not 0 < min_identity <= 1:
        raise MicrobiomeSuiteError(
            "--min-identity must be greater than 0 and less than or equal to 1."
        )
    table = ensure_input(table)
    rep_seqs = ensure_input(rep_seqs)
    database = _ensure_tax4fun2_database(database, database_mode=database_mode)
    _validate_tax4fun2_inputs(table, rep_seqs)
    resolved_threads = resolve_threads(threads)
    if output_dir.exists():
        if not output_dir.is_dir():
            raise MicrobiomeSuiteError(f"Output path exists and is not a directory: {output_dir}")
        if any(output_dir.iterdir()) and not force:
            raise MicrobiomeSuiteError(
                f"Output directory exists, pass --force to overwrite: {output_dir}"
            )
    output_dir.parent.mkdir(parents=True, exist_ok=True)

    params = {
        "tax4fun2_version": TAX4FUN2_VERSION,
        "database_mode": database_mode,
        "min_identity": min_identity,
        "normalize_by_copy_number": True,
        "normalize_pathways": normalize_pathways,
        "threads": resolved_threads,
    }
    outputs = {
        "functions": str(output_dir / TAX4FUN2_FUNCTIONS),
        "pathways": str(output_dir / TAX4FUN2_PATHWAYS),
        "coverage": str(output_dir / TAX4FUN2_COVERAGE),
        "manifest": str(output_dir / TAX4FUN2_MANIFEST),
    }

    with TemporaryDirectory(
        dir=output_dir.parent, prefix=".microsuite-tax4fun2-"
    ) as stage_temp_dir:
        staged_output = Path(stage_temp_dir) / "result"
        invoke_r_script(
            backend="tax4fun2",
            script_package="microsuite.functional.r",
            script_name="tax4fun2",
            resolve_image=resolve_functional_image,
            positional=[
                rep_seqs,
                table,
                database,
                str(resolved_threads),
                database_mode,
                str(min_identity),
                str(normalize_pathways).upper(),
                staged_output,
            ],
            runtime=runtime,
            image=image,
            engine=engine,
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(
                task="functional_profile",
                backend="tax4fun2",
                inputs={
                    "table": str(table),
                    "rep_seqs": str(rep_seqs),
                    "database": str(database),
                },
                outputs=outputs,
                params=params,
            ),
            local_missing_message=(
                "Tax4Fun2 requires external Rscript, Tax4Fun2 1.1.5, jsonlite, BLAST+, "
                "and compatible Tax4Fun2 v2 reference data. Install them locally, or use "
                "--runtime docker with the r-functional-tax4fun2 image."
            ),
        )
        _validate_tax4fun2_outputs(staged_output)
        sidecar = Path(stage_temp_dir) / "tax4fun2_container.json"
        if sidecar.exists():
            sidecar.replace(staged_output / sidecar.name)
        _replace_output_dir(staged_output, output_dir, force=force)


def _read_fasta_ids(
    path: Path, *, strict_headers: bool = False, reject_gzip: bool = False
) -> list[str]:
    if reject_gzip and path.name.lower().endswith(".gz"):
        raise MicrobiomeSuiteError(
            "PICRUSt2 representative FASTA must be plain text; gzipped FASTA is not supported "
            "by the upstream pipeline."
        )
    identifiers: list[str] = []
    seen: set[str] = set()
    has_sequence = False
    current: str | None = None
    allowed = set("ACGTRYSWKMBDHVN")
    if path.name.lower().endswith(".gz"):
        handle_context = gzip.open(path, "rt", encoding="utf-8")
    else:
        handle_context = path.open(encoding="utf-8")
    with handle_context as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None and not has_sequence:
                    raise MicrobiomeSuiteError(
                        f"Representative FASTA record has no sequence: {current}"
                    )
                header = line[1:].strip()
                if not header:
                    raise MicrobiomeSuiteError(
                        f"Representative FASTA has an empty identifier at line {line_number}."
                    )
                if strict_headers and any(character.isspace() for character in header):
                    raise MicrobiomeSuiteError(
                        "Representative FASTA headers must contain exactly one field without "
                        f"whitespace (line {line_number})."
                    )
                current = header if strict_headers else header.split(maxsplit=1)[0]
                if current in seen:
                    raise MicrobiomeSuiteError(
                        f"Representative FASTA identifier is duplicated: {current}"
                    )
                seen.add(current)
                identifiers.append(current)
                has_sequence = False
                continue
            if current is None:
                raise MicrobiomeSuiteError(
                    "Representative FASTA must begin with a '>' identifier line."
                )
            sequence = line.upper()
            invalid = sorted(set(sequence) - allowed)
            if invalid:
                raise MicrobiomeSuiteError(
                    f"Representative FASTA record {current} contains invalid DNA symbols: "
                    f"{''.join(invalid)}"
                )
            has_sequence = True
    if not identifiers:
        raise MicrobiomeSuiteError("Representative FASTA contains no records.")
    if not has_sequence:
        raise MicrobiomeSuiteError(f"Representative FASTA record has no sequence: {current}")
    return identifiers


def _read_tax4fun2_table_ids(path: Path, *, label: str = "Tax4Fun2") -> list[str]:
    with _open_table_text(path) as handle:
        rows = csv.reader(handle, delimiter="\t")
        try:
            header = next(rows)
        except StopIteration as exc:
            raise MicrobiomeSuiteError(f"{label} table is empty.") from exc
        if _is_mothur_shared_header(header):
            return _read_mothur_shared_ids(rows, header, label=label)
        if len(header) < 2:
            raise MicrobiomeSuiteError(
                f"{label} table must be tab-delimited with feature IDs in the first column "
                "and at least one sample column."
            )
        if not header[0].strip():
            raise MicrobiomeSuiteError(f"{label} table feature-ID column name cannot be empty.")
        sample_names = [name.strip() for name in header[1:]]
        if any(not name for name in sample_names) or len(sample_names) != len(set(sample_names)):
            raise MicrobiomeSuiteError(f"{label} sample names must be non-empty and unique.")

        identifiers: list[str] = []
        seen: set[str] = set()
        totals = [0.0] * len(sample_names)
        for line_number, row in enumerate(rows, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != len(header):
                raise MicrobiomeSuiteError(
                    f"{label} table line {line_number} has {len(row)} fields; "
                    f"expected {len(header)}."
                )
            identifier = row[0].strip()
            if not identifier:
                raise MicrobiomeSuiteError(
                    f"{label} table has an empty feature ID at line {line_number}."
                )
            if identifier in seen:
                raise MicrobiomeSuiteError(f"{label} feature ID is duplicated: {identifier}")
            seen.add(identifier)
            identifiers.append(identifier)
            for index, value in enumerate(row[1:]):
                try:
                    abundance = float(value)
                except ValueError as exc:
                    raise MicrobiomeSuiteError(
                        f"{label} abundance is not numeric at line {line_number}, "
                        f"sample {sample_names[index]}: {value!r}"
                    ) from exc
                if not math.isfinite(abundance) or abundance < 0:
                    raise MicrobiomeSuiteError(
                        f"{label} abundances must be finite and non-negative; "
                        f"found {value!r} for feature {identifier}, sample {sample_names[index]}."
                    )
                totals[index] += abundance
    if not identifiers:
        raise MicrobiomeSuiteError(f"{label} table contains no feature rows.")
    empty_samples = [name for name, total in zip(sample_names, totals, strict=True) if total <= 0]
    if empty_samples:
        raise MicrobiomeSuiteError(
            f"{label} samples must have positive total abundance: " + ", ".join(empty_samples)
        )
    return identifiers


def _open_table_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8-sig", newline="")
    return path.open(newline="", encoding="utf-8-sig")


def _is_mothur_shared_header(header: list[str]) -> bool:
    return len(header) >= 4 and [cell.strip().lower() for cell in header[:3]] == [
        "label",
        "group",
        "numotus",
    ]


def _read_mothur_shared_ids(rows: Any, header: list[str], *, label: str) -> list[str]:
    feature_ids = [cell.strip() for cell in header[3:]]
    if not feature_ids or any(not value for value in feature_ids):
        raise MicrobiomeSuiteError(f"{label} mothur shared table has empty feature IDs.")
    if len(feature_ids) != len(set(feature_ids)):
        raise MicrobiomeSuiteError(f"{label} mothur shared table feature IDs must be unique.")
    sample_names: list[str] = []
    totals: list[float] = []
    data_rows = 0
    for line_number, row in enumerate(rows, start=2):
        if not row or not any(cell.strip() for cell in row):
            continue
        if len(row) != len(header):
            raise MicrobiomeSuiteError(
                f"{label} mothur shared table line {line_number} has {len(row)} fields; "
                f"expected {len(header)}."
            )
        sample = row[1].strip()
        if not sample or sample in sample_names:
            raise MicrobiomeSuiteError(f"{label} mothur sample names must be non-empty and unique.")
        try:
            num_otus = int(row[2])
        except ValueError as exc:
            raise MicrobiomeSuiteError(
                f"{label} mothur shared table numOtus is not an integer at line {line_number}."
            ) from exc
        if num_otus != len(feature_ids):
            raise MicrobiomeSuiteError(
                f"{label} mothur shared table numOtus does not match its feature columns "
                f"at line {line_number}."
            )
        total = 0.0
        for value in row[3:]:
            try:
                abundance = float(value)
            except ValueError as exc:
                raise MicrobiomeSuiteError(
                    f"{label} abundance is not numeric at line {line_number}: {value!r}"
                ) from exc
            if not math.isfinite(abundance) or abundance < 0:
                raise MicrobiomeSuiteError(f"{label} abundances must be finite and non-negative.")
            total += abundance
        sample_names.append(sample)
        totals.append(total)
        data_rows += 1
    if not data_rows:
        raise MicrobiomeSuiteError(f"{label} mothur shared table contains no sample rows.")
    empty = [name for name, total in zip(sample_names, totals, strict=True) if total <= 0]
    if empty:
        raise MicrobiomeSuiteError(
            f"{label} samples must have positive total abundance: " + ", ".join(empty)
        )
    return feature_ids


def _read_biom_ids(path: Path, *, label: str) -> list[str]:
    try:
        biom_module = importlib.import_module("biom")
    except ImportError as exc:
        raise MicrobiomeSuiteError(
            f"{label} BIOM input requires the optional dependency 'biom-format'. "
            "Install with: uv sync --extra biom"
        ) from exc
    try:
        biom_table = biom_module.load_table(str(path))
        feature_ids = [str(value) for value in biom_table.ids(axis="observation")]
        sample_ids = [str(value) for value in biom_table.ids(axis="sample")]
        matrix = biom_table.matrix_data
    except Exception as exc:
        raise MicrobiomeSuiteError(f"{label} BIOM table is malformed: {path}") from exc
    if (
        not feature_ids
        or len(feature_ids) != len(set(feature_ids))
        or any(not value.strip() for value in feature_ids)
    ):
        raise MicrobiomeSuiteError(f"{label} BIOM feature IDs must be non-empty and unique.")
    if (
        not sample_ids
        or len(sample_ids) != len(set(sample_ids))
        or any(not value.strip() for value in sample_ids)
    ):
        raise MicrobiomeSuiteError(f"{label} BIOM sample IDs must be non-empty and unique.")
    try:
        totals = [0.0] * len(sample_ids)
        shape = getattr(matrix, "shape", None)
        expected_shape = (len(feature_ids), len(sample_ids))
        if shape is not None and tuple(shape) != expected_shape:
            raise MicrobiomeSuiteError(f"{label} BIOM matrix dimensions do not match its IDs.")
        sparse_matrix = matrix.tocsr() if hasattr(matrix, "tocsr") else None
        if sparse_matrix is not None:
            for row_index in range(len(feature_ids)):
                row = sparse_matrix.getrow(row_index)
                for index, value in zip(row.indices, row.data, strict=True):
                    abundance = float(value)
                    if not math.isfinite(abundance) or abundance < 0:
                        raise MicrobiomeSuiteError(
                            f"{label} abundances must be finite and non-negative."
                        )
                    totals[int(index)] += abundance
        else:
            for row in matrix:
                for index, value in enumerate(row):
                    abundance = float(value)
                    if not math.isfinite(abundance) or abundance < 0:
                        raise MicrobiomeSuiteError(
                            f"{label} abundances must be finite and non-negative."
                        )
                    totals[index] += abundance
    except MicrobiomeSuiteError:
        raise
    except (AttributeError, IndexError, TypeError, ValueError, OverflowError) as exc:
        raise MicrobiomeSuiteError(f"{label} BIOM table is malformed: {path}") from exc
    empty = [name for name, total in zip(sample_ids, totals, strict=True) if total <= 0]
    if empty:
        raise MicrobiomeSuiteError(
            f"{label} samples must have positive total abundance: " + ", ".join(empty)
        )
    return feature_ids


def _read_picrust2_table_ids(path: Path) -> list[str]:
    if path.suffix.lower() in {".biom", ".hdf5", ".h5"}:
        return _read_biom_ids(path, label="PICRUSt2")
    return _read_tax4fun2_table_ids(path, label="PICRUSt2")


def _validate_picrust2_inputs(table: Path, rep_seqs: Path) -> None:
    table_ids = _read_picrust2_table_ids(table)
    fasta_ids = _read_fasta_ids(rep_seqs, strict_headers=True, reject_gzip=True)
    if set(table_ids) != set(fasta_ids):
        missing_fasta = sorted(set(table_ids) - set(fasta_ids))
        missing_table = sorted(set(fasta_ids) - set(table_ids))
        details: list[str] = []
        if missing_fasta:
            details.append("missing from FASTA: " + ", ".join(missing_fasta[:5]))
        if missing_table:
            details.append("missing from table: " + ", ".join(missing_table[:5]))
        raise MicrobiomeSuiteError(
            "PICRUSt2 table and representative FASTA feature IDs must match exactly ("
            + "; ".join(details)
            + ")."
        )


def _validate_tax4fun2_inputs(table: Path, rep_seqs: Path) -> None:
    table_ids = _read_tax4fun2_table_ids(table)
    fasta_ids = _read_fasta_ids(rep_seqs)
    if set(table_ids) != set(fasta_ids):
        missing_fasta = sorted(set(table_ids) - set(fasta_ids))
        missing_table = sorted(set(fasta_ids) - set(table_ids))
        details: list[str] = []
        if missing_fasta:
            details.append("missing from FASTA: " + ", ".join(missing_fasta[:5]))
        if missing_table:
            details.append("missing from table: " + ", ".join(missing_table[:5]))
        raise MicrobiomeSuiteError(
            "Tax4Fun2 table and representative FASTA feature IDs must match exactly ("
            + "; ".join(details)
            + ")."
        )


def _ensure_tax4fun2_database(path: Path, *, database_mode: str) -> Path:
    path = _ensure_dir(path, flag="--database")
    required = (
        Path(database_mode) / f"{database_mode}.fasta",
        Path("KEGG") / "ko.txt",
        Path("KEGG") / "ko2ptw.txt",
        Path("KEGG") / "ptw.txt",
    )
    missing = [str(relative) for relative in required if not (path / relative).is_file()]
    profiles = path / database_mode
    if profiles.is_dir() and not any(profiles.glob("*.tbl.gz")):
        missing.append(f"{database_mode}/*.tbl.gz")
    if missing:
        raise MicrobiomeSuiteError(
            "Tax4Fun2 reference data is incomplete for "
            f"{database_mode}; missing: {', '.join(missing)}"
        )
    return path


def _validate_tax4fun2_outputs(output: Path) -> None:
    for filename in TAX4FUN2_REQUIRED_OUTPUTS:
        result = output / filename
        if not result.is_file() or result.stat().st_size == 0:
            raise MicrobiomeSuiteError(f"Tax4Fun2 did not produce its required output: {filename}")
    expected_headers = {
        TAX4FUN2_FUNCTIONS: "KO",
        TAX4FUN2_PATHWAYS: "pathway",
        TAX4FUN2_COVERAGE: "sample",
    }
    for filename, first_column in expected_headers.items():
        with (output / filename).open(encoding="utf-8") as handle:
            header = handle.readline().rstrip("\r\n").split("\t")
        if not header or header[0] != first_column:
            raise MicrobiomeSuiteError(
                f"Tax4Fun2 output {filename} has an invalid schema; "
                f"expected first column '{first_column}'."
            )


def _replace_output_dir(staged: Path, output: Path, *, force: bool) -> None:
    backup: Path | None = None
    if output.exists():
        if output.is_dir() and any(output.iterdir()) and not force:
            raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {output}")
        if not output.is_dir() and not force:
            raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {output}")
        backup = output.with_name(f".{output.name}.microsuite-backup-{uuid.uuid4().hex}")
        output.replace(backup)
    try:
        staged.replace(output)
    except OSError as exc:
        if output.exists():
            _remove_path(output)
        if backup is not None and backup.exists():
            backup.replace(output)
        raise MicrobiomeSuiteError(f"Failed to replace output directory: {output}") from exc
    if backup is not None and backup.exists():
        _remove_path(backup)


def _remove_path(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink()


def functional_profile_humann(
    *,
    reads: Path | None,
    database: Path | None,
    protein_database: Path | None,
    output_dir: Path,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if reads is None:
        raise MicrobiomeSuiteError("--reads is required for --backend humann.")
    executable = shutil.which("humann")
    if executable is None:
        raise MicrobiomeSuiteError(
            "HUMAnN functional profiling requires the external 'humann' command."
        )

    output_dir = prepare_output_dir(output_dir, force=force)
    command = [
        executable,
        "--input",
        str(ensure_input(reads)),
        "--output",
        str(output_dir),
        "--threads",
        str(resolve_threads(threads)),
    ]
    inputs = {"reads": str(reads)}
    if database is not None:
        command.extend(["--nucleotide-database", str(_ensure_dir(database, flag="--database"))])
        inputs["database"] = str(database)
    if protein_database is not None:
        command.extend(
            ["--protein-database", str(_ensure_dir(protein_database, flag="--protein-database"))]
        )
        inputs["protein_database"] = str(protein_database)

    run_command(
        command,
        "HUMAnN functional profiling failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="functional_profile",
            backend="humann",
            inputs=inputs,
            outputs={"output_dir": str(output_dir)},
        ),
    )


def prepare_output_dir(path: Path, *, force: bool) -> Path:
    if path.exists():
        if not path.is_dir():
            raise MicrobiomeSuiteError(f"Output path exists and is not a directory: {path}")
        if any(path.iterdir()) and not force:
            raise MicrobiomeSuiteError(
                f"Output directory exists, pass --force to overwrite: {path}"
            )
    path.mkdir(parents=True, exist_ok=True)
    return path


def _ensure_dir(path: Path, *, flag: str) -> Path:
    if not path.exists():
        raise MicrobiomeSuiteError(f"{flag} does not exist: {path}")
    if not path.is_dir():
        raise MicrobiomeSuiteError(f"{flag} must be a directory: {path}")
    return path
