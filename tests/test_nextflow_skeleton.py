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


def test_nextflow_docs_state_profiles_and_static_status() -> None:
    docs = (ROOT / "docs" / "api-nextflow.md").read_text(encoding="utf-8")

    assert "-profile local" in docs
    assert "-profile docker" in docs
    assert "-profile singularity" in docs
    assert "module files are placeholders" in docs
    assert "default tests validate the skeleton statically" in docs
