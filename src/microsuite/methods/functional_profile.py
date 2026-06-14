from __future__ import annotations

import shutil
from importlib.resources import files
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_BACKENDS = ("picrust2", "tax4fun2", "humann")

TAX4FUN2_SCRIPT = files("microsuite.functional.r").joinpath("tax4fun2.R")


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
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "functional_profile")
    if backend == "picrust2":
        functional_profile_picrust2(
            table=table,
            rep_seqs=rep_seqs,
            output_dir=output_dir,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
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
) -> None:
    if table is None:
        raise MicrobiomeSuiteError("--table is required for --backend picrust2.")
    if rep_seqs is None:
        raise MicrobiomeSuiteError("--rep-seqs is required for --backend picrust2.")
    executable = shutil.which("picrust2_pipeline.py")
    if executable is None:
        raise MicrobiomeSuiteError(
            "PICRUSt2 functional profiling requires 'picrust2_pipeline.py'. "
            "Install PICRUSt2 or run in a PICRUSt2 environment."
        )

    output_dir = prepare_output_dir(output_dir, force=force)
    command = [
        executable,
        "-s",
        str(ensure_input(rep_seqs)),
        "-i",
        str(ensure_input(table)),
        "-o",
        str(output_dir),
        "-p",
        str(resolve_threads(threads)),
    ]
    run_command(
        command,
        "PICRUSt2 functional profiling failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="functional_profile",
            backend="picrust2",
            inputs={"table": str(table), "rep_seqs": str(rep_seqs)},
            outputs={"output_dir": str(output_dir)},
        ),
    )


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
    rscript = shutil.which("Rscript")
    if rscript is None:
        raise MicrobiomeSuiteError(
            "Tax4Fun2 functional profiling requires external 'Rscript' with the Tax4Fun2 package."
        )

    output_dir = prepare_output_dir(output_dir, force=force)
    command = [
        rscript,
        str(TAX4FUN2_SCRIPT),
        str(ensure_input(rep_seqs)),
        str(ensure_input(table)),
        str(_ensure_dir(database, flag="--database")),
        str(output_dir),
        str(resolve_threads(threads)),
        database_mode,
        str(min_identity),
        str(normalize_pathways).upper(),
    ]
    run_command(
        command,
        "Tax4Fun2 functional profiling failed.",
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
            outputs={"output_dir": str(output_dir)},
        ),
    )


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
