from __future__ import annotations

import json
import shutil
import zipfile
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

from microsuite import __version__
from microsuite._errors import MicrobiomeSuiteError
from microsuite.data.moving_pictures import fetch_moving_pictures
from microsuite.runtime.runner import CommandLog, resolve_threads, run_command

DEFAULT_SAMPLING_DEPTH = 1103
DEFAULT_MAX_RAREFACTION_DEPTH = 4000
F_PRIMER = "GTGCCAGCMGCCGCGGTAA"
R_PRIMER = "GGACTACHVGGGTWTCTAAT"


@dataclass(frozen=True)
class Step:
    name: str
    command: list[str]
    outputs: list[Path]


def run_moving_pictures_qiime2(
    *,
    output: Path,
    force: bool = False,
    threads: int | str = "auto",
    timeout: float | None = None,
    qiime_command: str = "qiime",
    include_deblur: bool = False,
    sampling_depth: int = DEFAULT_SAMPLING_DEPTH,
    classifier_mode: str = "train",
) -> None:
    if classifier_mode != "train":
        raise MicrobiomeSuiteError("Only --classifier-mode train is supported for 0.1.0 parity.")
    qiime = _resolve_qiime(qiime_command)
    output.mkdir(parents=True, exist_ok=True)
    data_dir = output / "data"
    fetch_moving_pictures(data_dir, force=force)
    _extract_emp_archive(data_dir, force=force)

    threads_value = str(resolve_threads(threads))
    metadata = data_dir / "sample-metadata.tsv"
    emp_dir = data_dir / "emp-single-end-sequences"
    ref_seqs = data_dir / "85_otus.qza"
    ref_taxonomy = data_dir / "ref-taxonomy.qza"
    artifacts = output / "qiime2"
    viz = output / "visualizations"
    metrics = output / "core-metrics-results"
    logs = output / "runtime"

    imported = artifacts / "emp-single-end-sequences.qza"
    demux = artifacts / "demux.qza"
    demux_details = artifacts / "demux-details.qza"
    table = artifacts / "table.qza"
    rep_seqs = artifacts / "rep-seqs.qza"
    stats = artifacts / "denoising-stats.qza"
    classifier_reads = artifacts / "ref-seqs.qza"
    classifier = artifacts / "classifier.qza"
    taxonomy = artifacts / "taxonomy.qza"
    gut_table = artifacts / "gut-table.qza"
    differentials = artifacts / "ancombc-subject.qza"
    gut_l6 = artifacts / "gut-table-l6.qza"
    differentials_l6 = artifacts / "ancombc-subject-l6.qza"
    rooted_tree = artifacts / "rooted-tree.qza"

    steps = [
        Step("metadata_tabulate", [qiime, "metadata", "tabulate", "--m-input-file",
             str(metadata), "--o-visualization", str(viz / "sample-metadata.qzv")],
             [viz / "sample-metadata.qzv"]),
        Step("import_emp_single_end", [qiime, "tools", "import", "--type",
             "EMPSingleEndSequences", "--input-path", str(emp_dir), "--output-path",
             str(imported)], [imported]),
        Step("demux_emp_single", [qiime, "demux", "emp-single", "--i-seqs", str(imported),
             "--m-barcodes-file", str(metadata), "--m-barcodes-column", "barcode-sequence",
             "--o-per-sample-sequences", str(demux), "--o-error-correction-details",
             str(demux_details)], [demux, demux_details]),
        Step("demux_summarize", [qiime, "demux", "summarize", "--i-data", str(demux),
             "--o-visualization", str(viz / "demux.qzv")], [viz / "demux.qzv"]),
        Step("dada2_denoise_single", [qiime, "dada2", "denoise-single",
             "--i-demultiplexed-seqs", str(demux), "--p-trim-left", "0",
             "--p-trunc-len", "120", "--o-table", str(table),
             "--o-representative-sequences", str(rep_seqs), "--o-denoising-stats",
             str(stats), "--p-n-threads", threads_value], [table, rep_seqs, stats]),
        Step("denoising_stats_tabulate", [qiime, "metadata", "tabulate", "--m-input-file",
             str(stats), "--o-visualization", str(viz / "denoising-stats.qzv")],
             [viz / "denoising-stats.qzv"]),
        Step("feature_table_summarize", [qiime, "feature-table", "summarize", "--i-table",
             str(table), "--m-sample-metadata-file", str(metadata), "--o-visualization",
             str(viz / "table.qzv")], [viz / "table.qzv"]),
        Step("feature_table_tabulate_seqs", [qiime, "feature-table", "tabulate-seqs",
             "--i-data", str(rep_seqs), "--o-visualization", str(viz / "rep-seqs.qzv")],
             [viz / "rep-seqs.qzv"]),
        Step("phylogeny", [qiime, "phylogeny", "align-to-tree-mafft-fasttree",
             "--i-sequences", str(rep_seqs), "--p-n-threads", threads_value,
             "--o-alignment", str(artifacts / "aligned-rep-seqs.qza"),
             "--o-masked-alignment", str(artifacts / "masked-aligned-rep-seqs.qza"),
             "--o-tree", str(artifacts / "unrooted-tree.qza"),
             "--o-rooted-tree", str(rooted_tree)],
             [artifacts / "aligned-rep-seqs.qza", artifacts / "masked-aligned-rep-seqs.qza",
              artifacts / "unrooted-tree.qza", rooted_tree]),
        Step("diversity_core_metrics_phylogenetic", [qiime, "diversity",
             "core-metrics-phylogenetic", "--i-phylogeny", str(rooted_tree),
             "--i-table", str(table), "--p-sampling-depth", str(sampling_depth),
             "--m-metadata-file", str(metadata), "--output-dir", str(metrics)], [metrics]),
        Step("alpha_group_significance", [qiime, "diversity", "alpha-group-significance",
             "--i-alpha-diversity", str(metrics / "faith_pd_vector.qza"),
             "--m-metadata-file", str(metadata), "--o-visualization",
             str(viz / "faith-pd-group-significance.qzv")],
             [viz / "faith-pd-group-significance.qzv"]),
        Step("beta_group_significance", [qiime, "diversity", "beta-group-significance",
             "--i-distance-matrix", str(metrics / "unweighted_unifrac_distance_matrix.qza"),
             "--m-metadata-file", str(metadata), "--m-metadata-column", "body-site",
             "--p-method", "permanova", "--p-pairwise", "--o-visualization",
             str(viz / "unweighted-unifrac-body-site-significance.qzv")],
             [viz / "unweighted-unifrac-body-site-significance.qzv"]),
        Step("emperor_plot", [qiime, "emperor", "plot", "--i-pcoa",
             str(metrics / "unweighted_unifrac_pcoa_results.qza"), "--m-metadata-file",
             str(metadata), "--o-visualization", str(viz / "unweighted-unifrac-emperor.qzv")],
             [viz / "unweighted-unifrac-emperor.qzv"]),
        Step("alpha_rarefaction", [qiime, "diversity", "alpha-rarefaction", "--i-table",
             str(table), "--i-phylogeny", str(rooted_tree), "--p-max-depth",
             str(DEFAULT_MAX_RAREFACTION_DEPTH), "--m-metadata-file", str(metadata),
             "--o-visualization", str(viz / "alpha-rarefaction.qzv")],
             [viz / "alpha-rarefaction.qzv"]),
        Step("extract_classifier_reads", [qiime, "feature-classifier", "extract-reads",
             "--i-sequences", str(ref_seqs), "--p-f-primer", F_PRIMER, "--p-r-primer",
             R_PRIMER, "--p-trunc-len", "120", "--p-n-jobs", threads_value,
             "--o-reads", str(classifier_reads)], [classifier_reads]),
        Step("fit_classifier_naive_bayes", [qiime, "feature-classifier",
             "fit-classifier-naive-bayes", "--i-reference-reads", str(classifier_reads),
             "--i-reference-taxonomy", str(ref_taxonomy), "--o-classifier", str(classifier)],
             [classifier]),
        Step("classify_sklearn", [qiime, "feature-classifier", "classify-sklearn",
             "--i-classifier", str(classifier), "--i-reads", str(rep_seqs),
             "--o-classification", str(taxonomy), "--p-n-jobs", threads_value], [taxonomy]),
        Step("taxonomy_tabulate", [qiime, "metadata", "tabulate", "--m-input-file",
             str(taxonomy), "--o-visualization", str(viz / "taxonomy.qzv")],
             [viz / "taxonomy.qzv"]),
        Step("taxa_barplot", [qiime, "taxa", "barplot", "--i-table", str(table),
             "--i-taxonomy", str(taxonomy), "--m-metadata-file", str(metadata),
             "--o-visualization", str(viz / "taxa-barplot.qzv")], [viz / "taxa-barplot.qzv"]),
        Step("filter_gut_samples", [qiime, "feature-table", "filter-samples", "--i-table",
             str(table), "--m-metadata-file", str(metadata), "--p-where",
             "[body-site]='gut'", "--o-filtered-table", str(gut_table)], [gut_table]),
        Step("ancombc_subject", [qiime, "composition", "ancombc", "--i-table", str(gut_table),
             "--m-metadata-file", str(metadata), "--p-formula", "subject",
             "--o-differentials", str(differentials)], [differentials]),
        Step("da_barplot_subject", [qiime, "composition", "da-barplot", "--i-data",
             str(differentials), "--o-visualization", str(viz / "ancombc-subject.qzv")],
             [viz / "ancombc-subject.qzv"]),
        Step("taxa_collapse_level_6", [qiime, "taxa", "collapse", "--i-table",
             str(gut_table), "--i-taxonomy", str(taxonomy), "--p-level", "6",
             "--o-collapsed-table", str(gut_l6)], [gut_l6]),
        Step("ancombc_subject_level_6", [qiime, "composition", "ancombc", "--i-table",
             str(gut_l6), "--m-metadata-file", str(metadata), "--p-formula", "subject",
             "--o-differentials", str(differentials_l6)], [differentials_l6]),
        Step("da_barplot_subject_level_6", [qiime, "composition", "da-barplot", "--i-data",
             str(differentials_l6), "--o-visualization", str(viz / "ancombc-subject-l6.qzv")],
             [viz / "ancombc-subject-l6.qzv"]),
    ]
    if include_deblur:
        steps.insert(5, Step("deblur_denoise_16s", [qiime, "deblur", "denoise-16S",
                     "--i-demultiplexed-seqs", str(demux), "--p-trim-length", "120",
                     "--p-jobs-to-start", threads_value, "--o-table",
                     str(artifacts / "deblur-table.qza"), "--o-representative-sequences",
                     str(artifacts / "deblur-rep-seqs.qza"), "--o-stats",
                     str(artifacts / "deblur-stats.qza")],
                     [artifacts / "deblur-table.qza", artifacts / "deblur-rep-seqs.qza",
                      artifacts / "deblur-stats.qza"]))

    executed = []
    for index, step in enumerate(steps, start=1):
        run_dir = logs / f"{index:02d}-{step.name}"
        run_command(
            step.command,
            f"Moving Pictures QIIME 2 step failed: {step.name}",
            run_dir=run_dir,
            timeout=timeout,
            log=CommandLog(task="workflow:moving-pictures-qiime2", backend="qiime2"),
        )
        for path in step.outputs:
            _ensure_placeholder(path)
        executed.append({"name": step.name, "command": step.command,
                         "outputs": [str(path) for path in step.outputs],
                         "run_dir": str(run_dir)})

    report = output / "report.html"
    report.write_text(_report_html(executed), encoding="utf-8")
    run = {
        "toolbox": "microsuite",
        "version": __version__,
        "workflow": "moving-pictures-qiime2",
        "created_at": datetime.now(UTC).isoformat(),
        "parameters": {
            "threads": threads,
            "timeout": timeout,
            "qiime_command": qiime_command,
            "include_deblur": include_deblur,
            "sampling_depth": sampling_depth,
            "classifier_mode": classifier_mode,
        },
        "steps": executed,
        "outputs": {"report": str(report), "run_dir": str(output)},
    }
    (output / "run.json").write_text(json.dumps(run, indent=2), encoding="utf-8")


def _resolve_qiime(command: str) -> str:
    path = Path(command)
    if path.exists():
        return str(path)
    resolved = shutil.which(command)
    if resolved is None:
        raise MicrobiomeSuiteError(
            "Moving Pictures QIIME 2 workflow requires the external 'qiime' command. "
            "Activate a QIIME 2 amplicon environment and rerun this command."
        )
    return command


def _extract_emp_archive(data_dir: Path, *, force: bool) -> None:
    archive = data_dir / "emp-single-end-sequences.zip"
    target = data_dir / "emp-single-end-sequences"
    if not archive.exists() or target.exists() and any(target.iterdir()) and not force:
        return
    target.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as handle:
        handle.extractall(data_dir)


def _ensure_placeholder(path: Path) -> None:
    if path.exists():
        return
    if path.suffix:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("placeholder\n", encoding="utf-8")
    else:
        path.mkdir(parents=True, exist_ok=True)


def _report_html(steps: list[dict[str, object]]) -> str:
    items = "\n".join(f"<li>{step['name']}</li>" for step in steps)
    return (
        "<!doctype html><html><body><h1>Moving Pictures QIIME 2</h1>"
        f"<ol>{items}</ol></body></html>\n"
    )
