from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CONTAINERS = ROOT / "containers"


def test_container_dockerfiles_exist_with_expected_tools() -> None:
    expected = {
        "microsuite": ["microsuite", "uv"],
        "fastqc": ["fastqc", "openjdk"],
        "multiqc": ["multiqc"],
        "fastp": ["fastp"],
        "cutadapt": ["cutadapt"],
        "trimmomatic": ["trimmomatic", "openjdk"],
        "trim-galore": ["trim_galore", "cutadapt"],
        "vsearch": ["vsearch"],
        "usearch": ["usearch"],
        "mafft-fasttree": ["mafft", "FastTree"],
        "metaphlan": ["metaphlan"],
        "qiime2-amplicon": ["qiime", "QIIME 2"],
        "r-dada2": ["Rscript", "dada2"],
        "r-diffab-ancombc": ["Rscript", "ANCOMBC"],
        "r-diffab-aldex2": ["Rscript", "ALDEx2"],
        "r-diffab-maaslin2": ["Rscript", "Maaslin2"],
        "r-diffab-lefse": ["Rscript", "lefser"],
        "kraken2": ["kraken2", "bracken"],
        "microsuite-picrust2": ["microsuite", "picrust2"],
        "microsuite-dada2": ["microsuite", "Rscript", "dada2"],
    }

    for name, tokens in expected.items():
        dockerfile = CONTAINERS / name / "Dockerfile"
        assert dockerfile.exists(), name
        text = dockerfile.read_text(encoding="utf-8")
        assert "org.opencontainers.image.title" in text
        assert "org.opencontainers.image.description" in text
        assert "# Expected commands:" in text
        for token in tokens:
            assert token in text


def test_r_dada2_uses_digest_pinned_prebuilt_package() -> None:
    text = (CONTAINERS / "r-dada2" / "Dockerfile").read_text(encoding="utf-8")
    assert "quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2@sha256:" in text
    assert "packageVersion('dada2')) == '1.34.0'" in text
    assert "BiocManager::install" not in text


def test_microsuite_dada2_uses_digest_pinned_prebuilt_package() -> None:
    text = (CONTAINERS / "microsuite-dada2" / "Dockerfile").read_text(encoding="utf-8")
    assert "quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2@sha256:" in text
    assert "AS dada2-runtime" in text
    assert "FROM debian:bookworm-slim" in text
    assert "RUN rm -rf /usr/local/*" in text
    assert "COPY --from=dada2-runtime /usr/local/ /usr/local/" in text
    assert "packageVersion('dada2')) == '1.34.0'" in text
    assert "BiocManager::install" not in text


def test_microsuite_cli_images_support_arbitrary_runtime_users() -> None:
    base = (CONTAINERS / "microsuite" / "Dockerfile").read_text(encoding="utf-8")
    dada2 = (CONTAINERS / "microsuite-dada2" / "Dockerfile").read_text(encoding="utf-8")
    assert "chmod -R a+rX /opt/microsuite" in base
    assert "UV_PYTHON_INSTALL_DIR=/opt/uv/python" in dada2
    assert "chmod -R a+rX /opt/microsuite /opt/uv" in dada2


def test_nextflow_facing_images_install_process_metrics_tool() -> None:
    for image in ("cutadapt", "fastp", "microsuite", "microsuite-dada2", "multiqc"):
        text = (CONTAINERS / image / "Dockerfile").read_text(encoding="utf-8")
        assert "procps" in text, f"{image} must provide ps for Nextflow task metrics"


def test_microsuite_image_supports_downloaded_blast_binaries() -> None:
    text = (CONTAINERS / "microsuite" / "Dockerfile").read_text(encoding="utf-8")
    assert "libgomp1" in text


def test_r_diffab_smoke_outputs_use_writable_tmp() -> None:
    for backend in ("ancombc", "aldex2", "maaslin2", "lefse"):
        text = (CONTAINERS / f"r-diffab-{backend}" / "Dockerfile").read_text(encoding="utf-8")
        assert "/tmp/microsuite-smoke-out.tsv" in text
        assert "/opt/microsuite/smoke_out.tsv" not in text


def test_ancombc_backend_uses_cross_version_tse_input() -> None:
    source = (ROOT / "src/microsuite/diffab/r/ancombc.R").read_text(encoding="utf-8")
    assert "TreeSummarizedExperiment::TreeSummarizedExperiment" in source
    assert "meta_data = metadata" not in source


def test_dockerignore_excludes_local_artifacts() -> None:
    text = (ROOT / ".dockerignore").read_text(encoding="utf-8")

    for pattern in [".venv", ".git", "dist", "/runs", "/data", "/results", ".tokensave"]:
        assert pattern in text
    assert "src/microsuite/data" not in text


def test_methods_reference_links_backends_to_environments() -> None:
    text = (ROOT / "docs" / "methods.md").read_text(encoding="utf-8")

    assert "### Quality Reports" in text
    assert "### Quality Filtering" in text
    assert "### Feature Table Generation" in text
    assert "### Denoising And Clustering Backends" in text
    assert "### Differential Abundance" in text
    assert "| Subtopic |" not in text
    assert "Python function name" not in text
    assert "microsuite.api." not in text
    assert "microsuite.methods." not in text
    assert "microsuite.viz." not in text
    assert (
        "| Backend | Version | Status | CLI command | Python invocation | "
        "Image / environment | Operational tradeoff | Purpose |"
    ) in text
    assert "[microsuite Python](../containers/microsuite/Dockerfile)" in text
    assert "[FastQC](../containers/fastqc/Dockerfile)" in text
    assert "[MultiQC](../containers/multiqc/Dockerfile)" in text
    assert "[QIIME 2 amplicon](../containers/qiime2-amplicon/Dockerfile)" in text
    assert "[R diffab ancombc](../containers/r-diffab-ancombc/Dockerfile)" in text
    assert "[R diffab aldex2](../containers/r-diffab-aldex2/Dockerfile)" in text
    assert "[R diffab maaslin2](../containers/r-diffab-maaslin2/Dockerfile)" in text
    assert "[R diffab lefse](../containers/r-diffab-lefse/Dockerfile)" in text
    assert "[Kraken2](../containers/kraken2/Dockerfile)" in text
    assert "[VSEARCH](../containers/vsearch/Dockerfile)" in text
    assert "[USEARCH 12](../containers/usearch/Dockerfile)" in text
    assert "[MAFFT/FastTree](../containers/mafft-fasttree/Dockerfile)" in text
    assert "[MetaPhlAn](../containers/metaphlan/Dockerfile)" in text
    assert "[R DADA2](../containers/r-dada2/Dockerfile)" in text
    assert "| Generate ASV table | `qiime2-dada2` |" in text
    assert "| Generate ASV table | `qiime2-deblur` |" in text
    assert "| Generate ASV table | `dada2-r` |" in text
    assert "| Generate OTU-style table | `vsearch` |" in text
    assert "| Generate OTU-style table | `usearch` |" in text
    assert "| Generate QIIME-clustered table | `qiime2-vsearch` |" in text
    assert "ASV tables contain inferred exact sequence variants" in text
    assert "OTU-style tables contain" in text
    assert "features clustered by sequence identity" in text
    assert (
        "| `fastqc` | FastQC 0.12.1 | ready | "
        "`microsuite qc --backend fastqc` | "
        '`qc(backend="fastqc", inputs=[...], output_dir=...)` |'
    ) in text
    assert (
        "| `multiqc` | MultiQC user env | ready | "
        "`microsuite qc --backend multiqc` | "
        '`qc(backend="multiqc", input_dir=..., output_dir=...)` |'
    ) in text
    assert (
        "| `qiime2-demux-summarize` | QIIME 2 2024.10 | ready | "
        "`microsuite qc --backend qiime2-demux` | "
        '`qc(backend="qiime2-demux", demux=..., output=...)` |'
    ) in text
    assert (
        "| `qiime2-exclude-seqs` | QIIME 2 user env | ready | "
        "`microsuite qc_filter --backend qiime2-exclude-seqs` | "
        '`qc_filter(backend="qiime2-exclude-seqs", query_sequences=..., '
        "reference_sequences=...)` |"
    ) in text
    assert (
        "| `qiime2-bowtie2-build` | QIIME 2 user env | ready | "
        "`microsuite qc_filter --backend qiime2-bowtie2-build` | "
        '`qc_filter(backend="qiime2-bowtie2-build", sequences=..., output=...)` |'
    ) in text
    assert (
        "| `qiime2-decontam` | QIIME 2 user env | ready | "
        "`microsuite decontam --backend qiime2-decontam` | "
        '`decontam(backend="qiime2-decontam", table=..., metadata=..., output=...)` |'
    ) in text
    assert (
        "| `fastp` | fastp user env | ready | "
        "`microsuite trim --backend fastp` | "
        '`trim(backend="fastp", read1=..., output1=...)` |'
    ) in text
    assert (
        "| `cutadapt` | Cutadapt >=4.x user env | ready | "
        "`microsuite trim --backend cutadapt` | "
        '`trim(backend="cutadapt", read1=..., output1=..., adapter=...)` |'
    ) in text
    assert (
        "| `trimmomatic` | Trimmomatic >=0.39 user env | ready | "
        "`microsuite trim --backend trimmomatic` | "
        '`trim(backend="trimmomatic", read1=..., output1=..., trimmomatic_steps=[...])` |'
    ) in text
    assert (
        "| `trim-galore` | Trim Galore 0.6.x or v2.x user env | ready | "
        "`microsuite trim --backend trim-galore` | "
        '`trim(backend="trim-galore", read1=..., output1=..., trim_galore_version="auto")` |'
    ) in text
    assert (
        "| `qiime2-cutadapt` | QIIME 2 2024.10 | ready | "
        "`microsuite trim --backend qiime2-cutadapt --read1 demux.qza "
        "--output1 trimmed.qza` | "
        '`trim(backend="qiime2-cutadapt", read1=demux, output1=trimmed, adapter=...)` |'
    ) in text
    assert (
        "| `dada2-r` | DADA2 1.40.0 / Bioconductor 3.23 target | ready | "
        "`microsuite denoise --backend dada2-r` | "
        '`denoise(backend="dada2-r", demux=reads_dir, output_table=...)` |'
    ) in text
    assert (
        "| `kraken2` | Kraken2 2.1.3 | ready | "
        "`microsuite tax_classify --backend kraken2 --classifier DB` | "
        '`tax_classify(backend="kraken2", rep_seqs=..., classifier=database, output=...)` |'
    ) in text
    assert (
        "| `bracken` | Bracken user env | ready | "
        "`microsuite tax_classify --backend bracken --classifier DB --level S --read-length 150` | "
        '`tax_classify(backend="bracken", rep_seqs=kraken_report, '
        'classifier=database, output=..., level="S", read_length=150)` |'
    ) in text
    assert (
        "| `metaphlan` | MetaPhlAn user env | ready | "
        "`microsuite tax_classify --backend metaphlan --input-type fastq` | "
        '`tax_classify(backend="metaphlan", rep_seqs=..., output=..., input_type="fastq")` |'
    ) in text
    assert (
        "| `emu` | EMU user env | ready | "
        "`microsuite tax_classify --backend emu --input-type map-ont "
        "--classifier EMU_DB -o sample_rel-abundance.tsv` | "
        '`tax_classify(backend="emu", rep_seqs=long_reads, classifier=emu_db, '
        'output=..., input_type="map-ont")` |'
    ) in text
    assert (
        "| `mafft-fasttree` | MAFFT + FastTree user env | ready | "
        "`microsuite phylogeny --backend mafft-fasttree` | "
        '`phylogeny(backend="mafft-fasttree", rep_seqs=..., '
        "output_aligned=..., output_tree=...)` |"
    ) in text
    assert (
        "| `ancombc` | R image; ANCOMBC via Bioconductor | ready | "
        "`microsuite diff_abundance --backend ancombc` | "
        '`diff_abundance(backend="ancombc", table=..., group=..., output=...)` |'
    ) in text
    assert (
        "| `aldex2` | R image; ALDEx2 via Bioconductor | ready | "
        "`microsuite diff_abundance --backend aldex2` | "
        '`diff_abundance(backend="aldex2", table=..., group=..., output=...)` |'
    ) in text
    assert (
        "| `maaslin2` | R image; MaAsLin2 via Bioconductor | ready | "
        "`microsuite diff_abundance --backend maaslin2` | "
        '`diff_abundance(backend="maaslin2", table=..., group=..., output=...)` |'
    ) in text
    assert (
        "| `lefse` | R image; lefser via Bioconductor | ready | "
        "`microsuite diff_abundance --backend lefse` | "
        '`diff_abundance(backend="lefse", table=..., group=..., output=...)` |'
    ) in text
    assert (
        "| `native-beta-significance` | microsuite 0.1.0 | ready | "
        "`microsuite diversity beta-significance beta.tsv --metadata metadata.tsv "
        "--column body_site` | "
        '`beta_significance(distance_matrix, metadata, column=..., method="permanova")` |'
    ) in text
    assert (
        "| `mantel` | microsuite 0.1.0 | ready | "
        "`microsuite diversity mantel beta-a.tsv beta-b.tsv -o mantel.tsv` | "
        '`mantel_test(matrix_a, matrix_b, method="spearman")` |'
    ) in text
    assert (
        "| `rda` | microsuite 0.1.0 | ready | "
        "`microsuite diversity constrained-ordination table.h5ad --constraint body_site "
        "--method rda -o rda.tsv` | "
        '`constrained_ordination(adata, constraints=[...], method="rda")` |'
    ) in text
    assert (
        "| `native-gamma-diversity` | microsuite 0.1.0 | ready | "
        "`microsuite diversity gamma table.h5ad --group body_site -o gamma.tsv` | "
        '`gamma_diversity(adata, group=..., metric="observed_features")` |'
    ) in text
    assert (
        "| `beta-turnover` | microsuite 0.1.0 | ready | "
        "`microsuite diversity beta-turnover table.h5ad --level genus -o turnover.tsv` | "
        '`beta_turnover(adata, level="genus")` |'
    ) in text
    assert (
        "| `taxa-turnover` | microsuite 0.1.0 | ready | "
        "`microsuite diversity taxa-turnover table.h5ad --group body_site "
        "--level genus -o taxa-turnover.tsv` | "
        '`taxa_turnover(adata, group=..., level="genus")` |'
    ) in text
    assert (
        "| `native-correlation` | microsuite 0.1.0 | ready | "
        "`microsuite network infer --backend native-correlation --table table.h5ad "
        "-o network.tsv` | "
        '`network(backend="native-correlation", table=..., output=...)` |'
    ) in text
    assert (
        "| `sparcc` | microsuite 0.1.0 | ready | "
        "`microsuite network infer --backend sparcc --table table.h5ad -o sparcc.tsv` | "
        '`network(backend="sparcc", table=..., output=...)` |'
    ) in text
    assert (
        "| `spieceasi` | SpiecEasi R user env | ready | "
        "`microsuite network infer --backend spieceasi --table table.h5ad -o spieceasi.tsv` | "
        '`network(backend="spieceasi", table=..., output=...)` |'
    ) in text
    assert (
        "| `flashweave` | FlashWeave Julia user env | ready | "
        "`microsuite network infer --backend flashweave --table table.h5ad "
        "-o flashweave.edgelist` | "
        '`network(backend="flashweave", table=..., output=...)` |'
    ) in text
    assert (
        "| `randomforest` | microsuite 0.1.0 | ready | "
        "`microsuite ml classify --backend randomforest --table table.h5ad "
        "--target body_site -o predictions.tsv` | "
        '`ml_classify(backend="randomforest", table=..., target=..., output=...)` |'
    ) in text
    assert (
        "| `xgboost` | XGBoost optional Python package | ready | "
        "`microsuite ml classify --backend xgboost --table table.h5ad "
        "--target body_site -o predictions.tsv` | "
        '`ml_classify(backend="xgboost", table=..., target=..., output=...)` |'
    ) in text
    assert (
        "| `native-time-series` | microsuite 0.1.0 | ready | "
        "`microsuite ml longitudinal --backend native-time-series --table table.h5ad "
        "--subject subject --time day -o slopes.tsv` | "
        '`longitudinal(backend="native-time-series", table=..., subject=..., time=..., '
        "output=...)` |"
    ) in text
    assert (
        "| `picrust2` | PICRUSt2 user env | ready | "
        "`microsuite functional_profile --backend picrust2 --table table.biom "
        "--rep-seqs rep-seqs.fasta --output-dir functions` | "
        '`functional_profile(backend="picrust2", table=..., rep_seqs=..., output_dir=...)` |'
    ) in text
    assert (
        "| `tax4fun2` | Tax4Fun2 R user env | ready | "
        "`microsuite functional_profile --backend tax4fun2 --table otu-table.tsv "
        "--rep-seqs otus.fasta --database Tax4Fun2_ReferenceData_v2 --output-dir functions` | "
        '`functional_profile(backend="tax4fun2", table=..., rep_seqs=..., '
        "database=..., output_dir=...)` |"
    ) in text
    assert (
        "| `humann` | HUMAnN user env | ready | "
        "`microsuite functional_profile --backend humann --reads reads.fastq.gz "
        "--output-dir functions` | "
        '`functional_profile(backend="humann", reads=..., output_dir=...)` |'
    ) in text


@pytest.mark.skipif(
    os.environ.get("MICROSUITE_RUN_DOCKER_TESTS") != "1",
    reason="set MICROSUITE_RUN_DOCKER_TESTS=1 to run Docker-backed integration tests",
)
def test_fastqc_container_runs_tiny_fastq_end_to_end(tmp_path: Path) -> None:
    if shutil.which("docker") is None:
        pytest.skip("docker is not installed")

    image = "microsuite/fastqc:ci"
    image_check = subprocess.run(
        ["docker", "image", "inspect", image],
        check=False,
        text=True,
        capture_output=True,
        timeout=30,
    )
    if image_check.returncode != 0:
        pytest.skip(f"Docker image is not available locally: {image}")

    fixture = ROOT / "tests" / "fixtures" / "fastq" / "tiny.fastq"
    output_dir = tmp_path / "fastqc"
    output_dir.mkdir()

    command = [
        "docker",
        "run",
        "--rm",
        "-v",
        f"{fixture.parent.resolve()}:/input:ro",
        "-v",
        f"{output_dir.resolve()}:/output",
        image,
        "--outdir",
        "/output",
        "/input/tiny.fastq",
    ]
    try:
        result = subprocess.run(
            command,
            check=False,
            text=True,
            capture_output=True,
            timeout=90,
        )
    except subprocess.TimeoutExpired as exc:
        raise AssertionError(f"FastQC container timed out: {' '.join(command)}") from exc

    assert result.returncode == 0, result.stderr or result.stdout
    assert (output_dir / "tiny_fastqc.html").exists()
    assert (output_dir / "tiny_fastqc.zip").exists()
