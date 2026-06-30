from __future__ import annotations

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
NEXTFLOW = ROOT / "workflows" / "nextflow"


def test_nextflow_core_files_and_profiles_exist() -> None:
    expected = [
        NEXTFLOW / "main.nf",
        NEXTFLOW / "nextflow.config",
        NEXTFLOW / "profiles" / "local.config",
        NEXTFLOW / "profiles" / "docker.config",
        NEXTFLOW / "profiles" / "singularity.config",
    ]
    for path in expected:
        assert path.exists(), path


def test_nextflow_amplicon_modules_are_declared() -> None:
    modules = [
        "fastqc.nf",
        "multiqc.nf",
        "qiime2_dada2.nf",
        "qiime2_taxonomy.nf",
        "qiime2_phylogeny.nf",
        "qiime2_diversity.nf",
        "report.nf",
    ]
    for module in modules:
        path = NEXTFLOW / "modules" / module
        assert path.exists(), module

    main = (NEXTFLOW / "main.nf").read_text(encoding="utf-8")
    config = (NEXTFLOW / "nextflow.config").read_text(encoding="utf-8")
    assert "amplicon_qiime2" in main
    assert "manifest" in main
    assert "classifier" in main
    assert "MULTIQC.out.report_dir," not in main
    for module in modules:
        include_name = module.removesuffix(".nf")
        assert f"./modules/{include_name}" in main
    assert "profiles" in config
    for profile in ["local", "docker", "singularity"]:
        assert profile in config


def test_nextflow_amplicon_microsuite_workflow_drives_cli() -> None:
    modules = [
        "ms_cluster.nf",
        "ms_import.nf",
        "ms_diversity.nf",
        "ms_functional.nf",
        "ms_report.nf",
    ]
    for module in modules:
        assert (NEXTFLOW / "modules" / module).exists(), module

    main = (NEXTFLOW / "main.nf").read_text(encoding="utf-8")
    assert "amplicon_microsuite" in main
    assert "reads_fasta" in main
    for module in modules:
        assert f"./modules/{module.removesuffix('.nf')}" in main
    # every microsuite step is driven through the CLI
    for call in [
        "MS_CLUSTER",
        "MS_IMPORT",
        "MS_DIVERSITY",
        "MS_FUNCTIONAL",
        "MS_REPORT",
    ]:
        assert call in main


def test_nextflow_microsuite_modules_use_cli_and_stubs() -> None:
    expected_commands = {
        "ms_cluster.nf": ["microsuite cluster", "--reads", "vsearch --derep_fulllength", "stub:"],
        "ms_import.nf": ["microsuite import tsv", "stub:"],
        "ms_diversity.nf": ["microsuite diversity alpha", "breakaway", "inext", "stub:"],
        "ms_functional.nf": ["microsuite functional_profile", "picrust2", "stub:"],
        "ms_report.nf": ["report.html", "run.json", "stub:"],
    }
    for module, tokens in expected_commands.items():
        text = (NEXTFLOW / "modules" / module).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, f"{module}: {token}"
        assert "placeholder" not in text


def test_nextflow_modules_use_real_commands_and_stubs() -> None:
    expected_commands = {
        "fastqc.nf": ["fastqc --outdir", "stub:"],
        "multiqc.nf": ["multiqc", "stub:"],
        "qiime2_dada2.nf": ["qiime tools import", "qiime dada2 denoise", "stub:"],
        "qiime2_taxonomy.nf": ["qiime feature-classifier classify-sklearn", "stub:"],
        "qiime2_phylogeny.nf": ["qiime phylogeny align-to-tree-mafft-fasttree", "stub:"],
        "qiime2_diversity.nf": ["qiime diversity core-metrics-phylogenetic", "stub:"],
        "report.nf": ["report.html", "run.json", "stub:"],
    }

    for module, tokens in expected_commands.items():
        text = (NEXTFLOW / "modules" / module).read_text(encoding="utf-8")
        for token in tokens:
            assert token in text, f"{module}: {token}"
        assert "placeholder" not in text


def test_nextflow_profiles_assign_process_containers() -> None:
    docker = (NEXTFLOW / "profiles" / "docker.config").read_text(encoding="utf-8")
    singularity = (NEXTFLOW / "profiles" / "singularity.config").read_text(encoding="utf-8")

    for text in [docker, singularity]:
        for label in [
            "fastqc",
            "multiqc",
            "qiime2",
            "microsuite",
            "microsuite_amplicon",
            "microsuite_picrust2",
        ]:
            assert f"withLabel: {label}" in text
        assert "process.container =" not in text

    for image in [
        "microsuite/fastqc:ci",
        "microsuite/multiqc:ci",
        "microsuite/qiime2-amplicon:ci",
        "microsuite/microsuite:ci",
        "microsuite/prjna321534-alpha:ci",
        "microsuite/microsuite-picrust2:ci",
    ]:
        assert image in docker


def test_nextflow_docs_state_profiles_and_stub_status() -> None:
    docs = (ROOT / "docs" / "api-nextflow.md").read_text(encoding="utf-8")

    assert "-profile local" in docs
    assert "-profile docker" in docs
    assert "-profile singularity" in docs
    assert "module files contain runnable commands" in docs
    assert "Nextflow `-stub-run`" in docs
