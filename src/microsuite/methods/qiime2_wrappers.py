from __future__ import annotations

import json
import shutil
from datetime import UTC, datetime
from pathlib import Path

from microsuite import __version__
from microsuite._errors import MicrobiomeSuiteError
from microsuite._paths import ensure_input, prepare_output
from microsuite.methods._qiime import require_qiime, run_qiime
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

SUPPORTED_METHODS = {
    "metadata_tabulate": ("qiime2",),
    "qiime_import": ("qiime2-emp-single-end",),
    "demux": ("qiime2-emp-single",),
    "feature_summarize": ("qiime2",),
    "phylogeny": ("qiime2-mafft-fasttree", "mafft-fasttree"),
    "diversity_core": ("qiime2-core-metrics-phylogenetic",),
    "diversity_test": (
        "qiime2-alpha-group-significance",
        "qiime2-beta-group-significance",
        "qiime2-adonis",
    ),
    "ordination_plot": ("qiime2-emperor",),
    "rarefaction": ("qiime2-alpha-rarefaction",),
    "tax_train": ("qiime2-naive-bayes",),
    "tax_barplot": ("qiime2",),
    "feature_filter": ("qiime2-filter-samples",),
    "tax_collapse": ("qiime2",),
    "diff_viz": ("qiime2-da-barplot",),
}


def metadata_tabulate(
    *,
    backend: str,
    input_file: Path | None,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("metadata_tabulate", backend)
    input_file = _required(input_file, "--input-file", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 metadata tabulate")
    ensure_input(input_file)
    prepare_output(output, force=force)
    _run(
        [
            qiime,
            "metadata",
            "tabulate",
            "--m-input-file",
            str(input_file),
            "--o-visualization",
            str(output),
        ],
        "QIIME 2 metadata tabulate failed.",
        run_dir,
        timeout,
        "metadata_tabulate",
        backend,
    )


def qiime_import(
    *,
    backend: str,
    input_path: Path | None,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("qiime_import", backend)
    input_path = _required(input_path, "--input-path", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 EMP single-end import")
    if not input_path.exists():
        raise MicrobiomeSuiteError(f"Input path does not exist: {input_path}")
    prepare_output(output, force=force)
    _run(
        [
            qiime,
            "tools",
            "import",
            "--type",
            "EMPSingleEndSequences",
            "--input-path",
            str(input_path),
            "--output-path",
            str(output),
        ],
        "QIIME 2 EMP single-end import failed.",
        run_dir,
        timeout,
        "qiime_import",
        backend,
    )


def demux(
    *,
    backend: str,
    seqs: Path | None,
    metadata: Path | None,
    barcode_column: str | None,
    output_demux: Path | None,
    output_details: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("demux", backend)
    seqs = _required(seqs, "--seqs", backend)
    metadata = _required(metadata, "--metadata", backend)
    output_demux = _required(output_demux, "--output-demux", backend)
    output_details = _required(output_details, "--output-details", backend)
    if not barcode_column:
        raise MicrobiomeSuiteError(f"--barcode-column is required for --backend {backend}.")
    qiime = require_qiime("QIIME 2 EMP single-end demultiplexing")
    ensure_input(seqs)
    ensure_input(metadata)
    prepare_output(output_demux, force=force)
    prepare_output(output_details, force=force)
    _run(
        [
            qiime,
            "demux",
            "emp-single",
            "--i-seqs",
            str(seqs),
            "--m-barcodes-file",
            str(metadata),
            "--m-barcodes-column",
            barcode_column,
            "--o-per-sample-sequences",
            str(output_demux),
            "--o-error-correction-details",
            str(output_details),
        ],
        "QIIME 2 EMP single-end demultiplexing failed.",
        run_dir,
        timeout,
        "demux",
        backend,
    )


def feature_summarize(
    *,
    backend: str,
    mode: str,
    table: Path | None = None,
    rep_seqs: Path | None = None,
    metadata: Path | None = None,
    output: Path | None = None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("feature_summarize", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 feature-table visualization")
    prepare_output(output, force=force)
    if mode == "summarize":
        table = _required(table, "--table", backend)
        ensure_input(table)
        command = [qiime, "feature-table", "summarize", "--i-table", str(table)]
        if metadata is not None:
            ensure_input(metadata)
            command.extend(["--m-sample-metadata-file", str(metadata)])
        command.extend(["--o-visualization", str(output)])
    elif mode == "tabulate-seqs":
        rep_seqs = _required(rep_seqs, "--rep-seqs", backend)
        ensure_input(rep_seqs)
        command = [
            qiime,
            "feature-table",
            "tabulate-seqs",
            "--i-data",
            str(rep_seqs),
            "--o-visualization",
            str(output),
        ]
    else:
        raise MicrobiomeSuiteError("--mode must be 'summarize' or 'tabulate-seqs'.")
    _run(
        command,
        "QIIME 2 feature-table visualization failed.",
        run_dir,
        timeout,
        "feature_summarize",
        backend,
    )


def phylogeny(
    *,
    backend: str,
    rep_seqs: Path | None,
    output_aligned: Path | None,
    output_masked: Path | None,
    output_tree: Path | None,
    output_rooted_tree: Path | None,
    threads: int | str = 1,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    _expect_backend("phylogeny", backend)
    if backend == "mafft-fasttree":
        phylogeny_mafft_fasttree(
            rep_seqs=rep_seqs,
            output_aligned=output_aligned,
            output_masked=output_masked,
            output_tree=output_tree,
            output_rooted_tree=output_rooted_tree,
            threads=threads,
            force=force,
            run_dir=run_dir,
            timeout=timeout,
        )
        return
    rep_seqs = _required(rep_seqs, "--rep-seqs", backend)
    outputs = [
        _required(path, flag, backend)
        for path, flag in [
            (output_aligned, "--output-aligned"),
            (output_masked, "--output-masked"),
            (output_tree, "--output-tree"),
            (output_rooted_tree, "--output-rooted-tree"),
        ]
    ]
    qiime = require_qiime("QIIME 2 MAFFT/FastTree phylogeny")
    ensure_input(rep_seqs)
    for path in outputs:
        prepare_output(path, force=force)
    command = [
        qiime,
        "phylogeny",
        "align-to-tree-mafft-fasttree",
        "--i-sequences",
        str(rep_seqs),
        "--p-n-threads",
        str(resolve_threads(threads)),
        "--o-alignment",
        str(output_aligned),
        "--o-masked-alignment",
        str(output_masked),
        "--o-tree",
        str(output_tree),
        "--o-rooted-tree",
        str(output_rooted_tree),
    ]
    _run(command, "QIIME 2 phylogeny failed.", run_dir, timeout, "phylogeny", backend)


def phylogeny_mafft_fasttree(
    *,
    rep_seqs: Path | None,
    output_aligned: Path | None,
    output_masked: Path | None,
    output_tree: Path | None,
    output_rooted_tree: Path | None,
    threads: int | str,
    force: bool,
    run_dir: Path | None,
    timeout: float | None,
) -> None:
    rep_seqs = _required(rep_seqs, "--rep-seqs", "mafft-fasttree")
    output_aligned = _required(output_aligned, "--output-aligned", "mafft-fasttree")
    output_tree = _required(output_tree, "--output-tree", "mafft-fasttree")
    mafft = shutil.which("mafft")
    fasttree = shutil.which("FastTree") or shutil.which("fasttree")
    if mafft is None:
        raise MicrobiomeSuiteError(
            "Standalone phylogeny requires the external 'mafft' command. "
            "Install MAFFT or use the microsuite/mafft-fasttree container."
        )
    if fasttree is None:
        raise MicrobiomeSuiteError(
            "Standalone phylogeny requires the external 'FastTree' command. "
            "Install FastTree or use the microsuite/mafft-fasttree container."
        )

    ensure_input(rep_seqs)
    for path in (output_aligned, output_masked, output_tree, output_rooted_tree):
        if path is not None:
            prepare_output(path, force=force)

    mafft_command = [mafft, "--auto", "--thread", str(resolve_threads(threads)), str(rep_seqs)]
    mafft_result = run_command(
        mafft_command,
        "MAFFT alignment failed.",
        run_dir=run_dir / "mafft" if run_dir is not None else None,
        timeout=timeout,
        log=CommandLog(task="phylogeny", backend="mafft-fasttree"),
    )
    output_aligned.write_text(mafft_result.stdout, encoding="utf-8")
    if output_masked is not None:
        output_masked.write_text(mafft_result.stdout, encoding="utf-8")

    fasttree_command = [fasttree, "-nt", str(output_aligned)]
    fasttree_result = run_command(
        fasttree_command,
        "FastTree tree construction failed.",
        run_dir=run_dir / "fasttree" if run_dir is not None else None,
        timeout=timeout,
        log=CommandLog(task="phylogeny", backend="mafft-fasttree"),
    )
    output_tree.write_text(fasttree_result.stdout, encoding="utf-8")
    if output_rooted_tree is not None:
        output_rooted_tree.write_text(fasttree_result.stdout, encoding="utf-8")
    if run_dir is not None:
        _write_summary_run(
            run_dir,
            "phylogeny",
            "mafft-fasttree",
            [mafft_command, fasttree_command],
            {
                "aligned": output_aligned,
                "masked": output_masked or output_aligned,
                "tree": output_tree,
                "rooted_tree": output_rooted_tree or output_tree,
            },
        )


def diversity_core(
    *,
    backend: str,
    table: Path | None,
    phylogeny_path: Path | None,
    metadata: Path | None,
    sampling_depth: int,
    output_dir: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("diversity_core", backend)
    table = _required(table, "--table", backend)
    phylogeny_path = _required(phylogeny_path, "--phylogeny", backend)
    metadata = _required(metadata, "--metadata", backend)
    output_dir = _required(output_dir, "--output-dir", backend)
    if sampling_depth < 1:
        raise MicrobiomeSuiteError("--sampling-depth must be at least 1.")
    qiime = require_qiime("QIIME 2 core metrics phylogenetic")
    for path in (table, phylogeny_path, metadata):
        ensure_input(path)
    _prepare_directory_output(output_dir, force=force)
    _run(
        [
            qiime,
            "diversity",
            "core-metrics-phylogenetic",
            "--i-phylogeny",
            str(phylogeny_path),
            "--i-table",
            str(table),
            "--p-sampling-depth",
            str(sampling_depth),
            "--m-metadata-file",
            str(metadata),
            "--output-dir",
            str(output_dir),
        ],
        "QIIME 2 core metrics phylogenetic failed.",
        run_dir,
        timeout,
        "diversity_core",
        backend,
    )


def diversity_test(
    *,
    backend: str,
    alpha_diversity: Path | None = None,
    distance_matrix: Path | None = None,
    metadata: Path | None = None,
    metadata_column: str | None = None,
    output: Path | None = None,
    method: str = "permanova",
    pairwise: bool = False,
    formula: str | None = None,
    permutations: int = 999,
    n_jobs: int | str = 1,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    backend = backend.lower()
    _expect_backend("diversity_test", backend)
    metadata = _required(metadata, "--metadata", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 diversity group significance")
    ensure_input(metadata)
    prepare_output(output, force=force)
    if backend == "qiime2-adonis":
        distance_matrix = _required(distance_matrix, "--distance-matrix", backend)
        if not formula or not formula.strip():
            raise MicrobiomeSuiteError(f"--formula is required for --backend {backend}.")
        if "~" in formula or "\n" in formula or "\r" in formula:
            raise MicrobiomeSuiteError(
                "--formula must contain only the R formula right-hand side for QIIME 2."
            )
        if method.lower() != "permanova":
            raise MicrobiomeSuiteError(
                "--method must remain permanova for the QIIME 2 adonis backend."
            )
        if metadata_column:
            raise MicrobiomeSuiteError(
                f"--metadata-column is not used for --backend {backend}; use --formula."
            )
        if pairwise:
            raise MicrobiomeSuiteError(f"--pairwise is not supported for --backend {backend}.")
        if permutations < 1:
            raise MicrobiomeSuiteError("--permutations must be at least 1 for QIIME 2 adonis.")
        ensure_input(distance_matrix)
        command = [
            qiime,
            "diversity",
            "adonis",
            "--i-distance-matrix",
            str(distance_matrix),
            "--m-metadata-file",
            str(metadata),
            "--p-formula",
            formula.strip(),
            "--p-permutations",
            str(permutations),
            "--p-n-jobs",
            str(resolve_threads(n_jobs)),
            "--o-visualization",
            str(output),
        ]
    elif backend == "qiime2-alpha-group-significance":
        alpha_diversity = _required(alpha_diversity, "--alpha-diversity", backend)
        ensure_input(alpha_diversity)
        command = [
            qiime,
            "diversity",
            "alpha-group-significance",
            "--i-alpha-diversity",
            str(alpha_diversity),
            "--m-metadata-file",
            str(metadata),
            "--o-visualization",
            str(output),
        ]
    else:
        distance_matrix = _required(distance_matrix, "--distance-matrix", backend)
        if not metadata_column:
            raise MicrobiomeSuiteError(f"--metadata-column is required for --backend {backend}.")
        if formula:
            raise MicrobiomeSuiteError(
                f"--formula is only supported for --backend qiime2-adonis, not {backend}."
            )
        if permutations < 1:
            raise MicrobiomeSuiteError("--permutations must be at least 1 for QIIME 2.")
        ensure_input(distance_matrix)
        command = [
            qiime,
            "diversity",
            "beta-group-significance",
            "--i-distance-matrix",
            str(distance_matrix),
            "--m-metadata-file",
            str(metadata),
            "--m-metadata-column",
            metadata_column,
            "--p-method",
            method,
            "--p-permutations",
            str(permutations),
            "--o-visualization",
            str(output),
        ]
        command.append("--p-pairwise" if pairwise else "--p-no-pairwise")
    _run(
        command,
        "QIIME 2 diversity group significance failed.",
        run_dir,
        timeout,
        "diversity_test",
        backend,
    )


def ordination_plot(
    *,
    backend: str,
    pcoa: Path | None,
    metadata: Path | None,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("ordination_plot", backend)
    pcoa = _required(pcoa, "--pcoa", backend)
    metadata = _required(metadata, "--metadata", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 Emperor plot")
    ensure_input(pcoa)
    ensure_input(metadata)
    prepare_output(output, force=force)
    _run(
        [
            qiime,
            "emperor",
            "plot",
            "--i-pcoa",
            str(pcoa),
            "--m-metadata-file",
            str(metadata),
            "--o-visualization",
            str(output),
        ],
        "QIIME 2 Emperor plot failed.",
        run_dir,
        timeout,
        "ordination_plot",
        backend,
    )


def rarefaction(
    *,
    backend: str,
    table: Path | None,
    phylogeny_path: Path | None,
    metadata: Path | None,
    max_depth: int,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("rarefaction", backend)
    table = _required(table, "--table", backend)
    metadata = _required(metadata, "--metadata", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 alpha rarefaction")
    ensure_input(table)
    ensure_input(metadata)
    if phylogeny_path is not None:
        ensure_input(phylogeny_path)
    prepare_output(output, force=force)
    command = [
        qiime,
        "diversity",
        "alpha-rarefaction",
        "--i-table",
        str(table),
        "--p-max-depth",
        str(max_depth),
        "--m-metadata-file",
        str(metadata),
    ]
    if phylogeny_path is not None:
        command.extend(["--i-phylogeny", str(phylogeny_path)])
    command.extend(["--o-visualization", str(output)])
    _run(command, "QIIME 2 alpha rarefaction failed.", run_dir, timeout, "rarefaction", backend)


def tax_train(
    *,
    backend: str,
    ref_seqs: Path | None,
    ref_taxonomy: Path | None,
    f_primer: str | None,
    r_primer: str | None,
    trunc_len: int,
    output: Path | None,
    threads: int | str = 1,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("tax_train", backend)
    ref_seqs = _required(ref_seqs, "--ref-seqs", backend)
    ref_taxonomy = _required(ref_taxonomy, "--ref-taxonomy", backend)
    output = _required(output, "--output", backend)
    if not f_primer or not r_primer:
        raise MicrobiomeSuiteError("--f-primer and --r-primer are required.")
    qiime = require_qiime("QIIME 2 naive Bayes classifier training")
    ensure_input(ref_seqs)
    ensure_input(ref_taxonomy)
    prepare_output(output, force=force)
    reads = output.with_name(f"{output.stem}-reads.qza")
    prepare_output(reads, force=force)
    threads_value = str(resolve_threads(threads))
    first = [
        qiime,
        "feature-classifier",
        "extract-reads",
        "--i-sequences",
        str(ref_seqs),
        "--p-f-primer",
        f_primer,
        "--p-r-primer",
        r_primer,
        "--p-trunc-len",
        str(trunc_len),
        "--p-n-jobs",
        threads_value,
        "--o-reads",
        str(reads),
    ]
    second = [
        qiime,
        "feature-classifier",
        "fit-classifier-naive-bayes",
        "--i-reference-reads",
        str(reads),
        "--i-reference-taxonomy",
        str(ref_taxonomy),
        "--o-classifier",
        str(output),
    ]
    if run_dir is None:
        run_qiime(
            first,
            "QIIME 2 extract-reads failed.",
            timeout=timeout,
            task="tax_train",
            backend=backend,
        )
        run_qiime(
            second,
            "QIIME 2 fit-classifier-naive-bayes failed.",
            timeout=timeout,
            task="tax_train",
            backend=backend,
        )
    else:
        run_qiime(
            first,
            "QIIME 2 extract-reads failed.",
            run_dir=run_dir / "extract-reads",
            timeout=timeout,
            task="tax_train",
            backend=backend,
        )
        run_qiime(
            second,
            "QIIME 2 fit-classifier-naive-bayes failed.",
            run_dir=run_dir / "fit-classifier-naive-bayes",
            timeout=timeout,
            task="tax_train",
            backend=backend,
        )
        _write_summary_run(run_dir, "tax_train", backend, [first, second], {"classifier": output})


def tax_barplot(
    *,
    backend: str,
    table: Path | None,
    taxonomy: Path | None,
    metadata: Path | None,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("tax_barplot", backend)
    table = _required(table, "--table", backend)
    taxonomy = _required(taxonomy, "--taxonomy", backend)
    metadata = _required(metadata, "--metadata", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 taxa barplot")
    for path in (table, taxonomy, metadata):
        ensure_input(path)
    prepare_output(output, force=force)
    _run(
        [
            qiime,
            "taxa",
            "barplot",
            "--i-table",
            str(table),
            "--i-taxonomy",
            str(taxonomy),
            "--m-metadata-file",
            str(metadata),
            "--o-visualization",
            str(output),
        ],
        "QIIME 2 taxa barplot failed.",
        run_dir,
        timeout,
        "tax_barplot",
        backend,
    )


def feature_filter(
    *,
    backend: str,
    table: Path | None,
    metadata: Path | None,
    where: str | None,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("feature_filter", backend)
    table = _required(table, "--table", backend)
    metadata = _required(metadata, "--metadata", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 feature-table filter-samples")
    ensure_input(table)
    ensure_input(metadata)
    prepare_output(output, force=force)
    command = [
        qiime,
        "feature-table",
        "filter-samples",
        "--i-table",
        str(table),
        "--m-metadata-file",
        str(metadata),
        "--o-filtered-table",
        str(output),
    ]
    if where:
        command.extend(["--p-where", where])
    _run(command, "QIIME 2 filter-samples failed.", run_dir, timeout, "feature_filter", backend)


def tax_collapse(
    *,
    backend: str,
    table: Path | None,
    taxonomy: Path | None,
    level: int,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("tax_collapse", backend)
    table = _required(table, "--table", backend)
    taxonomy = _required(taxonomy, "--taxonomy", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 taxa collapse")
    ensure_input(table)
    ensure_input(taxonomy)
    prepare_output(output, force=force)
    _run(
        [
            qiime,
            "taxa",
            "collapse",
            "--i-table",
            str(table),
            "--i-taxonomy",
            str(taxonomy),
            "--p-level",
            str(level),
            "--o-collapsed-table",
            str(output),
        ],
        "QIIME 2 taxa collapse failed.",
        run_dir,
        timeout,
        "tax_collapse",
        backend,
    )


def diff_viz(
    *,
    backend: str,
    data: Path | None,
    output: Path | None,
    force: bool = False,
    run_dir: Path | None = None,
    timeout: float | None = None,
) -> None:
    _expect_backend("diff_viz", backend)
    data = _required(data, "--data", backend)
    output = _required(output, "--output", backend)
    qiime = require_qiime("QIIME 2 composition DA barplot")
    ensure_input(data)
    prepare_output(output, force=force)
    _run(
        [
            qiime,
            "composition",
            "da-barplot",
            "--i-data",
            str(data),
            "--o-visualization",
            str(output),
        ],
        "QIIME 2 DA barplot failed.",
        run_dir,
        timeout,
        "diff_viz",
        backend,
    )


def _expect_backend(method: str, backend: str) -> None:
    backend = backend.lower()
    if backend not in SUPPORTED_METHODS[method]:
        choices = ", ".join(SUPPORTED_METHODS[method])
        raise MicrobiomeSuiteError(
            f"Unsupported {method} backend '{backend}'. Choose one of: {choices}"
        )


def _required(value: Path | None, flag: str, backend: str) -> Path:
    if value is None:
        raise MicrobiomeSuiteError(f"{flag} is required for --backend {backend}.")
    return value


def _prepare_directory_output(path: Path, *, force: bool) -> None:
    if path.exists() and not path.is_dir():
        raise MicrobiomeSuiteError(f"Output path exists and is not a directory: {path}")
    if path.exists() and any(path.iterdir()) and not force:
        raise MicrobiomeSuiteError(f"Output directory is not empty, pass --force to use it: {path}")
    path.mkdir(parents=True, exist_ok=True)


def _run(
    command: list[str],
    failure_message: str,
    run_dir: Path | None,
    timeout: float | None,
    task: str,
    backend: str,
) -> None:
    run_qiime(
        command,
        failure_message,
        run_dir=run_dir,
        timeout=timeout,
        task=task,
        backend=backend.lower(),
    )


def _write_summary_run(
    run_dir: Path, task: str, backend: str, commands: list[list[str]], outputs: dict[str, Path]
) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    payload = {
        "toolbox": "microsuite",
        "version": __version__,
        "created_at": datetime.now(UTC).isoformat(),
        "task": task,
        "backend": backend,
        "commands": commands,
        "outputs": {key: str(value) for key, value in outputs.items()},
    }
    (run_dir / "run.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
