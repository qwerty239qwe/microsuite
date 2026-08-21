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
    assert "nf-core/setup-nextflow@v3.0.1" in text
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
    assert "mothur-smoke:" in text
    assert "mamba-org/setup-micromamba@v3.2.1" in text
    assert "mothur=1.48.5" in text
    assert 'MICROSUITE_RUN_MOTHUR_SMOKE: "1"' in text
    assert "uv run pytest tests/integration/test_mothur_smoke.py" in text
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
    assert "group: ${{ github.workflow }}-${{ github.ref }}-${{ github.event_name }}" in text
    assert "cancel-in-progress: true" in text
    assert "docker-build:" in text
    assert "docker/setup-buildx-action@v4.2.0" in text
    assert "docker/build-push-action@v7.3.0" in text
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
    assert "containers/r-diffab-ancombc/Dockerfile" in text
    assert "containers/r-diffab-aldex2/Dockerfile" in text
    assert "containers/r-diffab-maaslin2/Dockerfile" in text
    assert "containers/r-diffab-maaslin3/Dockerfile" in text
    assert "containers/r-diffab-lefse/Dockerfile" in text
    assert "containers/r-functional-tax4fun2/Dockerfile" in text
    assert "containers/r-ecology/Dockerfile" in text
    for basename in (
        "r-batch-mmuphin",
        "r-batch-combatseq",
        "r-batch-conqur",
        "r-batch-plsdabatch",
        "r-batch-metadict",
    ):
        assert f"containers/{basename}/Dockerfile" in text
    assert "containers/r-diffab/Dockerfile" not in text
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
    assert "image: r-diffab-ancombc" in text
    assert "image: r-diffab-aldex2" in text
    assert "image: r-diffab-maaslin2" in text
    assert "image: r-diffab-maaslin3" in text
    assert "image: r-diffab-lefse" in text
    assert "image: r-functional-tax4fun2" in text
    assert "image: r-batch-mmuphin" in text
    assert "image: r-batch-combatseq" in text
    assert "image: r-batch-conqur" in text
    assert "image: r-batch-plsdabatch" in text
    assert "containers/r-batch-metadict/Dockerfile" in text
    assert "image: r-batch-metadict" in text
    assert "image: microsuite-picrust2" in text
    assert "image: microsuite-dada2" in text
    assert "image: mothur" in text
    assert text.count("heavy: false") == 12
    assert text.count("heavy: true") == 17
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
    assert "Smoke test microsuite-dada2 image" in text
    assert "Smoke test PICRUSt2 image contents" in text
    assert "Provision swap for PICRUSt2-SC smoke" in text
    assert "sudo fallocate -l 16G /mnt/picrust2-swapfile" in text
    assert "sudo swapon /mnt/picrust2-swapfile" in text
    assert "Smoke test PICRUSt2 SC functional profiling" in text
    assert "Smoke test PICRUSt2 custom oldIMG single-reference functional profiling" in text
    assert text.index("Smoke test PICRUSt2 custom oldIMG") < text.index(
        "Provision swap for PICRUSt2-SC smoke"
    )
    assert "--entrypoint picrust2_pipeline.py microsuite/${{ matrix.image }}:ci --version" in text
    assert "--entrypoint python microsuite/${{ matrix.image }}:ci -c" in text
    assert "--picrust2-database SC" in text
    assert "--picrust2-coverage" in text
    assert "--entrypoint sh microsuite/${{ matrix.image }}:ci -lc" in text
    assert r"default_tables[\"EC\"]" in text
    assert r"default_tables[\"16S\"]" in text
    assert "/opt/conda/bin/python -c" in text
    assert "/opt/conda/bin/microsuite functional_profile" in text
    assert "--picrust2-database custom" in text
    assert "--picrust2-no-pathways" in text
    assert '--picrust2-ref-dir1 "$ref_dir"' in text
    assert '--picrust2-custom-trait-tables-ref1 "$ec_table"' in text
    assert '--picrust2-marker-gene-table-ref1 "$marker_table"' in text
    assert "custom_metagenome_count=$(find tmp/docker-picrust2-custom/picrust2" in text
    assert "-path '*/*_metagenome_out/pred_metagenome_unstrat.tsv.gz'" in text
    assert 'test "$custom_metagenome_count" -eq 1' in text
    assert "tmp/docker-picrust2-custom/picrust2/picrust2_manifest.json" in text
    assert "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz" in text
    assert "EC_metagenome_out/weighted_nsti.tsv.gz" in text
    assert "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" in text
    assert "pathways_out/path_abun_unstrat.tsv.gz" in text
    assert "pathways_out/path_cov_unstrat.tsv.gz" in text
    assert "picrust2_manifest.json" in text
    assert "containers/microsuite-dada2/Dockerfile" in text
    assert "denoise --backend dada2-r" in text
    assert "--user 10001:10001" in text
    assert "--entrypoint /opt/microsuite/.venv/bin/microsuite" in text
    assert "test -s tmp/docker-dada2/out/rep-seqs.fasta" in text
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
    assert "docker/login-action@v4.6.0" in text
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
