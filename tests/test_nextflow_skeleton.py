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
    for module in modules:
        include_name = module.removesuffix(".nf")
        assert f"./modules/{include_name}" in main
    assert "profiles" in config
    for profile in ["local", "docker", "singularity"]:
        assert profile in config


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


def test_nextflow_docs_state_profiles_and_stub_status() -> None:
    docs = (ROOT / "docs" / "api-nextflow.md").read_text(encoding="utf-8")

    assert "-profile local" in docs
    assert "-profile docker" in docs
    assert "-profile singularity" in docs
    assert "module files contain runnable commands" in docs
    assert "Nextflow `-stub-run`" in docs
