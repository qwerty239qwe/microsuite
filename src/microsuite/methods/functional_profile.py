from __future__ import annotations

import csv
import math
import shutil
from importlib.resources import files
from pathlib import Path
from tempfile import TemporaryDirectory

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input
from microsuite.methods._dispatch import require_backend
from microsuite.runtime.container import resolve_functional_image
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
) -> None:
    backend = require_backend(backend, SUPPORTED_BACKENDS, "functional_profile")
    if backend != "tax4fun2" and (runtime != "local" or image is not None or engine != "docker"):
        raise MicrobiomeSuiteError(
            "--runtime, --image, and --engine currently apply only to --backend tax4fun2."
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


def _read_fasta_ids(path: Path) -> list[str]:
    identifiers: list[str] = []
    seen: set[str] = set()
    has_sequence = False
    current: str | None = None
    allowed = set("ACGTRYSWKMBDHVN")
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None and not has_sequence:
                    raise MicrobiomeSuiteError(
                        f"Representative FASTA record has no sequence: {current}"
                    )
                current = line[1:].split(maxsplit=1)[0] if line[1:].strip() else ""
                if not current:
                    raise MicrobiomeSuiteError(
                        f"Representative FASTA has an empty identifier at line {line_number}."
                    )
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


def _read_tax4fun2_table_ids(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = csv.reader(handle, delimiter="\t")
        try:
            header = next(rows)
        except StopIteration as exc:
            raise MicrobiomeSuiteError("Tax4Fun2 table is empty.") from exc
        if len(header) < 2:
            raise MicrobiomeSuiteError(
                "Tax4Fun2 table must be tab-delimited with feature IDs in the first column "
                "and at least one sample column."
            )
        if not header[0].strip():
            raise MicrobiomeSuiteError("Tax4Fun2 table feature-ID column name cannot be empty.")
        sample_names = [name.strip() for name in header[1:]]
        if any(not name for name in sample_names) or len(sample_names) != len(set(sample_names)):
            raise MicrobiomeSuiteError("Tax4Fun2 sample names must be non-empty and unique.")

        identifiers: list[str] = []
        seen: set[str] = set()
        totals = [0.0] * len(sample_names)
        for line_number, row in enumerate(rows, start=2):
            if not row or not any(cell.strip() for cell in row):
                continue
            if len(row) != len(header):
                raise MicrobiomeSuiteError(
                    f"Tax4Fun2 table line {line_number} has {len(row)} fields; "
                    f"expected {len(header)}."
                )
            identifier = row[0].strip()
            if not identifier:
                raise MicrobiomeSuiteError(
                    f"Tax4Fun2 table has an empty feature ID at line {line_number}."
                )
            if identifier in seen:
                raise MicrobiomeSuiteError(f"Tax4Fun2 feature ID is duplicated: {identifier}")
            seen.add(identifier)
            identifiers.append(identifier)
            for index, value in enumerate(row[1:]):
                try:
                    abundance = float(value)
                except ValueError as exc:
                    raise MicrobiomeSuiteError(
                        f"Tax4Fun2 abundance is not numeric at line {line_number}, "
                        f"sample {sample_names[index]}: {value!r}"
                    ) from exc
                if not math.isfinite(abundance) or abundance < 0:
                    raise MicrobiomeSuiteError(
                        "Tax4Fun2 abundances must be finite and non-negative; "
                        f"found {value!r} for feature {identifier}, sample {sample_names[index]}."
                    )
                totals[index] += abundance
    if not identifiers:
        raise MicrobiomeSuiteError("Tax4Fun2 table contains no feature rows.")
    empty_samples = [name for name, total in zip(sample_names, totals, strict=True) if total <= 0]
    if empty_samples:
        raise MicrobiomeSuiteError(
            "Tax4Fun2 samples must have positive total abundance: " + ", ".join(empty_samples)
        )
    return identifiers


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
    if output.exists():
        if output.is_dir() and not any(output.iterdir()):
            output.rmdir()
        elif not force:
            raise MicrobiomeSuiteError(f"Output exists, pass --force to overwrite: {output}")
        elif output.is_dir():
            shutil.rmtree(output)
        else:
            output.unlink()
    staged.replace(output)


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
