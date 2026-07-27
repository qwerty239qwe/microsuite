from __future__ import annotations

import shutil
from pathlib import Path

from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.metadata import Artifact, stage_execution
from microsuite.methods._dispatch import reject_options, require_backend
from microsuite.methods.mothur import run_mothur, select_output
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_METHODS = ("qiime2", "kraken2", "bracken", "metaphlan", "emu", "mothur")


def tax_classify(
    *,
    backend: str,
    rep_seqs: Path,
    output: Path,
    classifier: Path | str | None = None,
    input_type: str = "fastq",
    level: str = "S",
    read_length: int = 150,
    taxonomy_reference: Path | None = None,
    taxonomy_map: Path | None = None,
    otu_list: Path | None = None,
    count_table: Path | None = None,
    threads: int | str = 1,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = require_backend(backend, SUPPORTED_METHODS, "taxonomy classification")

    if backend == "mothur":
        reject_options(
            "mothur",
            {
                "--classifier": classifier,
                # input_type defaults to "fastq" rather than None, so only a value the
                # user actually chose can be rejected here. mothur classifies FASTA and
                # has no input-type switch, so anything non-default is meaningless to it.
                "--input-type": None if input_type == "fastq" else input_type,
            },
        )
        if taxonomy_reference is None or taxonomy_map is None:
            raise MicrobiomeSuiteError(
                "--taxonomy-reference and --taxonomy-map are both required for --backend mothur."
            )
    else:
        reject_options(
            backend,
            {
                "--taxonomy-reference": taxonomy_reference,
                "--taxonomy-map": taxonomy_map,
                "--otu-list": otu_list,
                "--count-table": count_table,
            },
        )

    from microsuite.refdb.service import resolve_classifier

    if classifier is not None:
        classifier = resolve_classifier(str(classifier))

    resolved_threads = resolve_threads(threads)
    with stage_execution(
        output.resolve().parent,
        stage="taxonomy",
        backend=backend,
        params={"level": level, "threads": resolved_threads},
        outputs=[Artifact("taxonomy", output.resolve(), kind="taxonomy_table", required=False)],
    ):
        if backend == "qiime2":
            tax_classify_qiime2(
                rep_seqs=rep_seqs,
                classifier=classifier,
                output=output,
                threads=resolved_threads,
                force=force,
                run_dir=run_dir,
                timeout=timeout,
            )
            return
        if backend == "bracken":
            tax_classify_bracken(
                kraken_report=rep_seqs,
                database=classifier,
                output=output,
                level=level,
                read_length=read_length,
                force=force,
                run_dir=run_dir,
                timeout=timeout,
            )
            return
        if backend == "kraken2":
            tax_classify_kraken2(
                reads=rep_seqs,
                database=classifier,
                output=output,
                threads=resolved_threads,
                force=force,
                run_dir=run_dir,
                timeout=timeout,
            )
            return
        if backend == "metaphlan":
            tax_classify_metaphlan(
                reads=rep_seqs,
                database=classifier,
                output=output,
                input_type=input_type,
                threads=resolved_threads,
                force=force,
                run_dir=run_dir,
                timeout=timeout,
            )
            return
        if backend == "emu":
            tax_classify_emu(
                reads=rep_seqs,
                database=classifier,
                output=output,
                input_type=input_type,
                threads=resolved_threads,
                force=force,
                run_dir=run_dir,
                timeout=timeout,
            )
            return
        if backend == "mothur":
            # taxonomy_reference/taxonomy_map are validated non-None above, in
            # the `backend == "mothur"` branch before `stage_execution`. The
            # asserts restate that runtime guarantee for the type checker here
            # and, unlike a cast, fail loudly at this exact spot if a future
            # refactor ever breaks the dominance between validation and
            # dispatch, instead of passing None into ensure_input() silently.
            assert taxonomy_reference is not None
            assert taxonomy_map is not None
            tax_classify_mothur(
                rep_seqs=rep_seqs,
                output=output,
                taxonomy_reference=taxonomy_reference,
                taxonomy_map=taxonomy_map,
                otu_list=otu_list,
                count_table=count_table,
                force=force,
                run_dir=run_dir,
                timeout=timeout,
            )
            return


def tax_classify_emu(
    *,
    reads: Path,
    database: Path | None,
    output: Path,
    input_type: str,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    emu_type = "map-ont" if input_type == "fastq" else input_type
    if emu_type not in {"map-ont", "map-pb", "lr:hq", "map-hifi", "sr"}:
        raise MicrobiomeSuiteError(
            "--input-type must be one of: map-ont, map-pb, lr:hq, map-hifi, sr "
            "for --backend emu. The default CLI value 'fastq' maps to map-ont."
        )
    if not output.name.endswith("_rel-abundance.tsv"):
        raise MicrobiomeSuiteError(
            "--output must end with '_rel-abundance.tsv' for --backend emu so it matches "
            "the EMU output naming convention."
        )
    emu = shutil.which("emu")
    if emu is None:
        raise MicrobiomeSuiteError(
            "EMU long-read amplicon profiling requires the external 'emu' command. "
            "Install EMU and provide an EMU database with --classifier or EMU_DATABASE_DIR."
        )

    ensure_input(reads)
    if database is not None and not database.exists():
        raise MicrobiomeSuiteError(f"EMU database does not exist: {database}")
    prepare_output(output, force=force)

    basename = output.name[: -len("_rel-abundance.tsv")]
    command = [
        emu,
        "abundance",
        str(reads),
        "--type",
        emu_type,
        "--threads",
        str(threads),
        "--output-dir",
        str(output.parent),
        "--output-basename",
        basename,
    ]
    if database is not None:
        command.extend(["--db", str(database)])
    run_command(
        command,
        "EMU long-read amplicon profiling failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="tax_classify",
            backend="emu",
            inputs={"reads": str(reads), **({"database": str(database)} if database else {})},
            outputs={"relative_abundance": str(output)},
            params={"input_type": emu_type},
        ),
    )


def tax_classify_bracken(
    *,
    kraken_report: Path,
    database: Path | None,
    output: Path,
    level: str,
    read_length: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if database is None:
        raise MicrobiomeSuiteError("--classifier is required for --backend bracken.")
    if level not in {"D", "P", "C", "O", "F", "G", "S"}:
        raise MicrobiomeSuiteError(
            "--level must be one of: D, P, C, O, F, G, S for --backend bracken."
        )
    if read_length < 1:
        raise MicrobiomeSuiteError("--read-length must be at least 1 for --backend bracken.")
    bracken = shutil.which("bracken")
    if bracken is None:
        raise MicrobiomeSuiteError(
            "Bracken abundance re-estimation requires the external 'bracken' command. "
            "Install Bracken or use the microsuite/kraken2 container and rerun this command."
        )

    ensure_input(kraken_report)
    if not database.exists():
        raise MicrobiomeSuiteError(f"Bracken database does not exist: {database}")
    prepare_output(output, force=force)
    command = [
        bracken,
        "-d",
        str(database),
        "-i",
        str(kraken_report),
        "-o",
        str(output),
        "-r",
        str(read_length),
        "-l",
        level,
    ]
    run_command(
        command,
        "Bracken abundance re-estimation failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="tax_classify",
            backend="bracken",
            inputs={"kraken_report": str(kraken_report), "database": str(database)},
            outputs={"abundance": str(output)},
            params={"level": level, "read_length": read_length},
        ),
    )


def tax_classify_metaphlan(
    *,
    reads: Path,
    database: Path | None,
    output: Path,
    input_type: str,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if input_type not in {"fastq", "fasta", "bowtie2out", "sam"}:
        raise MicrobiomeSuiteError(
            "--input-type must be one of: fastq, fasta, bowtie2out, sam for --backend metaphlan."
        )
    metaphlan = shutil.which("metaphlan")
    if metaphlan is None:
        raise MicrobiomeSuiteError(
            "MetaPhlAn taxonomy profiling requires the external 'metaphlan' command. "
            "Install MetaPhlAn or use the microsuite/metaphlan container and rerun this command."
        )

    ensure_input(reads)
    if database is not None and not database.exists():
        raise MicrobiomeSuiteError(f"MetaPhlAn database does not exist: {database}")
    prepare_output(output, force=force)
    bowtie2out = output.with_suffix(".bowtie2.bz2")
    prepare_output(bowtie2out, force=force)

    command = [
        metaphlan,
        str(reads),
        "--input_type",
        input_type,
        "--nproc",
        str(threads),
        "--bowtie2out",
        str(bowtie2out),
        "-o",
        str(output),
    ]
    if database is not None:
        command.extend(["--bowtie2db", str(database)])
    run_command(
        command,
        "MetaPhlAn taxonomy profiling failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="tax_classify",
            backend="metaphlan",
            inputs={"reads": str(reads), **({"database": str(database)} if database else {})},
            outputs={"profile": str(output), "bowtie2out": str(bowtie2out)},
            params={"input_type": input_type},
        ),
    )


def tax_classify_kraken2(
    *,
    reads: Path,
    database: Path | None,
    output: Path,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if database is None:
        raise MicrobiomeSuiteError("--classifier is required for --backend kraken2.")
    kraken2 = shutil.which("kraken2")
    if kraken2 is None:
        raise MicrobiomeSuiteError(
            "Kraken2 taxonomy classification requires the external 'kraken2' command. "
            "Install Kraken2 or use the microsuite/kraken2 container and rerun this command."
        )

    ensure_input(reads)
    if not database.exists():
        raise MicrobiomeSuiteError(f"Kraken2 database does not exist: {database}")
    prepare_output(output, force=force)
    per_read_output = output.with_suffix(".kraken")
    prepare_output(per_read_output, force=force)

    command = [
        kraken2,
        "--db",
        str(database),
        "--threads",
        str(threads),
        "--report",
        str(output),
        "--output",
        str(per_read_output),
        str(reads),
    ]
    run_command(
        command,
        "Kraken2 taxonomy classification failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="tax_classify",
            backend="kraken2",
            inputs={"reads": str(reads), "database": str(database)},
            outputs={"report": str(output), "per_read": str(per_read_output)},
        ),
    )


def tax_classify_qiime2(
    *,
    rep_seqs: Path,
    classifier: Path | None,
    output: Path,
    threads: int,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    if classifier is None:
        raise MicrobiomeSuiteError("--classifier is required for --backend qiime2.")
    qiime = shutil.which("qiime")
    if qiime is None:
        raise MicrobiomeSuiteError(
            "QIIME 2 taxonomy classification requires the external 'qiime' command. "
            "Activate a QIIME 2 environment and rerun this command."
        )

    ensure_input(rep_seqs)
    ensure_input(classifier)
    prepare_output(output, force=force)

    command = [
        qiime,
        "feature-classifier",
        "classify-sklearn",
        "--i-classifier",
        str(classifier),
        "--i-reads",
        str(rep_seqs),
        "--o-classification",
        str(output),
        "--p-n-jobs",
        str(threads),
    ]
    run_command(
        command,
        "QIIME 2 classification failed.",
        run_dir=run_dir,
        timeout=timeout,
        log=CommandLog(
            task="tax_classify",
            backend="qiime2",
            inputs={"rep_seqs": str(rep_seqs), "classifier": str(classifier)},
            outputs={"taxonomy": str(output)},
        ),
    )


def tax_classify_mothur(
    *,
    rep_seqs: Path,
    output: Path,
    taxonomy_reference: Path,
    taxonomy_map: Path,
    otu_list: Path | None,
    count_table: Path | None,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    """Classify sequences with mothur's naive Bayes classifier.

    With an ``otu_list`` from the mothur cluster backend, the per-sequence
    assignments are consensus-collapsed per OTU via classify.otu.
    """
    ensure_input(rep_seqs)
    ensure_input(taxonomy_reference)
    ensure_input(taxonomy_map)
    prepare_output(output, force=force)

    work_dir = output.parent / f"{output.stem}.mothur"
    work_dir.mkdir(parents=True, exist_ok=True)

    def step(name: str, params: dict[str, str]) -> list[Path]:
        step_run_dir = None if run_dir is None else run_dir / name.replace(".", "_")
        return run_mothur(name, params, work_dir=work_dir, run_dir=step_run_dir, timeout=timeout)

    classify_params = {
        "fasta": str(rep_seqs),
        "reference": str(taxonomy_reference),
        "taxonomy": str(taxonomy_map),
    }
    if count_table is not None:
        ensure_input(count_table)
        classify_params["count"] = str(count_table)

    out = step("classify.seqs", classify_params)
    assignments = select_output(out, ".taxonomy", step="classify.seqs")

    if otu_list is not None:
        ensure_input(otu_list)
        otu_params = {"list": str(otu_list), "taxonomy": str(assignments)}
        if count_table is not None:
            otu_params["count"] = str(count_table)
        out = step("classify.otu", otu_params)
        assignments = select_output(out, ".cons.taxonomy", step="classify.otu")

    shutil.copyfile(assignments, output)
