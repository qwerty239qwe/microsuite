from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_github_actions_ci_runs_project_quality_gate() -> None:
    workflow = ROOT / ".github" / "workflows" / "ci.yml"
    assert workflow.exists()

    text = workflow.read_text(encoding="utf-8")
    assert "astral-sh/setup-uv" in text
    assert 'python-version: ["3.11", "3.12"]' in text
    assert "UV_PYTHON: ${{ matrix.python-version }}" in text
    assert "enable-cache: true" in text
    assert "cache-dependency-glob: uv.lock" in text
    assert "uv sync --all-extras --locked" in text
    assert "--extra dev --all-extras" not in text
    assert "uv run pytest" in text
    assert "uv run ruff check ." in text
    assert "uv run ruff format --check ." in text
    assert "uv run ty check" in text
    assert "uv build" in text
    assert "nextflow-smoke:" in text
    assert "nf-core/setup-nextflow@v2" in text
    assert "nextflow-io/setup-nextflow" not in text
    assert 'version: "25.10.0"' in text
    assert "nextflow run -stub-run workflows/nextflow/main.nf" in text
    assert "-stub-run" in text
    assert "--manifest tests/fixtures/nextflow/manifest.tsv" in text
    assert "--classifier tests/fixtures/nextflow/classifier.qza" in text
    assert "--outdir results/nextflow-smoke" in text
    assert "external-integration:" in text
    assert "External tool Python API integration" in text
    assert 'MICROSUITE_RUN_EXTERNAL_INTEGRATION: "1"' in text
    assert "sudo apt-get install -y --no-install-recommends fastp vsearch mafft fasttree" in text
    assert "uv pip install cutadapt multiqc" in text
    assert "uv run pytest tests/integration/test_external_tools.py" in text
    assert "docker-build:" not in text
    # A real (non-stub) end-to-end run is available on demand.
    assert "nextflow-real:" in text
    assert "run-real-nextflow" in text
    assert "nextflow run workflows/nextflow/main.nf" in text
    assert "results/nextflow-real/report/report.html" in text


def test_github_actions_docker_workflow_builds_and_tests_images() -> None:
    workflow = ROOT / ".github" / "workflows" / "docker.yml"
    assert workflow.exists()

    text = workflow.read_text(encoding="utf-8")
    assert "docker-build:" in text
    assert "docker/setup-buildx-action@v3" in text
    assert "docker/build-push-action@v6" in text
    assert "containers/microsuite/Dockerfile" in text
    assert "containers/fastqc/Dockerfile" in text
    assert "containers/multiqc/Dockerfile" in text
    assert "containers/fastp/Dockerfile" in text
    assert "containers/cutadapt/Dockerfile" in text
    assert "containers/trimmomatic/Dockerfile" in text
    assert "containers/trim-galore/Dockerfile" in text
    assert "containers/vsearch/Dockerfile" in text
    assert "containers/usearch/Dockerfile" in text
    assert "containers/mafft-fasttree/Dockerfile" in text
    assert "containers/metaphlan/Dockerfile" in text
    assert "containers/kraken2/Dockerfile" in text
    assert "containers/qiime2-amplicon/Dockerfile" in text
    assert "containers/r-dada2/Dockerfile" in text
    assert "containers/r-diffab/Dockerfile" in text
    assert "build-heavy-containers" in text
    assert "Skipping heavy container image" in text
    assert "image: microsuite" in text
    assert "image: fastqc" in text
    assert "image: multiqc" in text
    assert "image: fastp" in text
    assert "image: cutadapt" in text
    assert "image: trimmomatic" in text
    assert "image: trim-galore" in text
    assert "image: vsearch" in text
    assert "image: usearch" in text
    assert "image: mafft-fasttree" in text
    assert "image: metaphlan" in text
    assert "image: kraken2" in text
    assert "image: qiime2-amplicon" in text
    assert "image: r-dada2" in text
    assert "image: r-diffab" in text
    assert "image: microsuite-picrust2" in text
    assert "image: microsuite-dada2" in text
    assert text.count("heavy: false") == 11
    assert text.count("heavy: true") == 6
    assert "load: true" in text
    assert "Run FastQC tiny example" in text
    assert "Run MultiQC tiny example" in text
    assert "Run fastp tiny example" in text
    assert "Run Cutadapt tiny example" in text
    assert "Run Trimmomatic tiny example" in text
    assert "Run Trim Galore tiny example" in text
    assert "Run VSEARCH tiny example" in text
    assert "Run USEARCH tiny example" in text
    assert "Smoke test MAFFT/FastTree image" in text
    assert "Smoke test MetaPhlAn image" in text
    assert "Smoke test R diffab image" in text
    assert "Smoke test microsuite-dada2 image" in text
    assert "containers/microsuite-dada2/Dockerfile" in text
    assert "denoise --backend dada2-r" in text
    assert "test -s tmp/docker-dada2/rep-seqs.fasta" in text
    assert 'c("ANCOMBC", "ALDEx2", "Maaslin2", "lefser")' in text
    assert 'c("ancombc", "aldex2", "maaslin2", "lefse")' in text
    assert "--cluster_fast /input/tiny.fasta" in text
    assert "--minseqlength 1" in text
    assert "-cluster_fast /input/tiny.fasta" in text
    assert "test -s tmp/docker-fastp/tiny.fastq" in text
    assert "test -s tmp/docker-multiqc/output/multiqc_report.html" in text
    assert (
        "docker run --rm --entrypoint kraken2 microsuite/${{ matrix.image }}:ci --version" in text
    )
    assert "docker run --rm --entrypoint bracken microsuite/${{ matrix.image }}:ci -h" in text
    assert "command -v mafft && command -v FastTree" in text
    assert "microsuite/${{ matrix.image }}:ci --version" in text
    # Built + smoke-tested images are published to GHCR (main only).
    assert "packages: write" in text
    assert "docker/login-action@v3" in text
    assert "ghcr.io/${owner}/microsuite/${{ matrix.image }}" in text
    assert 'sha="sha-$(git rev-parse --short HEAD)"' in text
    assert "docker push" in text


def test_source_data_package_is_not_ignored() -> None:
    assert (ROOT / "src" / "microsuite" / "data" / "__init__.py").exists()
    assert (ROOT / "src" / "microsuite" / "data" / "moving_pictures.py").exists()
    assert (
        ROOT / "src" / "microsuite" / "data" / "fixtures" / "moving_pictures_small" / "table.tsv"
    ).exists()


def test_demo_data_attribution_is_documented() -> None:
    attribution = (ROOT / "docs" / "data-attribution.md").read_text(encoding="utf-8")
    fixture_readme = (
        ROOT / "src" / "microsuite" / "data" / "fixtures" / "moving_pictures_small" / "README.md"
    ).read_text(encoding="utf-8")

    for text in (attribution, fixture_readme):
        assert "not owned by this project" in text
        assert "10.1186/gb-2011-12-5-r50" in text
        assert "QIIME 2 Moving Pictures tutorial" in text
