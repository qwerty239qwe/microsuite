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
    assert "include { MS_FUNCTIONAL_OVERRIDES } from './modules/ms_functional'" in main
    assert "include { MS_FUNCTIONAL_SYMBOLIC_OVERRIDE } from './modules/ms_functional'" in main
    assert "include { MS_FUNCTIONAL_CUSTOM } from './modules/ms_functional'" in main
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


def test_nextflow_picrust2_defaults_and_safe_custom_arguments() -> None:
    module = (NEXTFLOW / "modules" / "ms_functional.nf").read_text(encoding="utf-8")
    config = (NEXTFLOW / "nextflow.config").read_text(encoding="utf-8")

    assert "'--picrust2-database', params.picrust2_database ?: 'SC'" in module
    for option in (
        "--picrust2-ref-dir1",
        "--picrust2-ref-dir2",
        "--picrust2-custom-trait-tables-ref1",
        "--picrust2-custom-trait-tables-ref2",
        "--picrust2-marker-gene-table-ref1",
        "--picrust2-marker-gene-table-ref2",
        "--picrust2-pathway-map",
        "--picrust2-reaction-func",
        "--picrust2-regroup-map",
        "--picrust2-max-nsti",
        "--picrust2-no-pathways",
        "--picrust2-coverage",
        "--picrust2-no-regroup",
    ):
        assert option in module
    assert "shell_quote" in module
    assert "picrust2_args" not in module
    assert "process MS_FUNCTIONAL_OVERRIDES" in module
    assert "path override_assets, stageAs: 'picrust2_override_asset??/*'" in module
    assert "val override_plan" in module
    assert "process MS_FUNCTIONAL_SYMBOLIC_OVERRIDE" in module
    assert "picrust2_database = 'SC'" in config
    for key in (
        "picrust2_ref_dir1",
        "picrust2_ref_dir2",
        "picrust2_custom_trait_tables_ref1",
        "picrust2_custom_trait_tables_ref2",
        "picrust2_marker_gene_table_ref1",
        "picrust2_marker_gene_table_ref2",
        "picrust2_pathway_map",
        "picrust2_reaction_func",
        "picrust2_regroup_map",
        "picrust2_max_nsti",
        "picrust2_no_pathways",
        "picrust2_coverage",
        "picrust2_no_regroup",
    ):
        assert key in config


def test_nextflow_custom_picrust2_assets_are_real_staged_path_inputs() -> None:
    module = (NEXTFLOW / "modules" / "ms_functional.nf").read_text(encoding="utf-8")
    main = (NEXTFLOW / "main.nf").read_text(encoding="utf-8")

    # The custom process receives one collection of path inputs. Nextflow
    # stages every member; custom_plan contains only indexes into that staged
    # collection, which prevents host paths from leaking into containers.
    assert "process MS_FUNCTIONAL_CUSTOM" in module
    assert "path custom_assets, stageAs: 'picrust2_custom_asset??/*'" in module
    assert "val custom_plan" in module
    assert "def staged_at" in module
    assert "custom_assets[index as int]" in module
    for option in (
        "staged_at(custom_plan.ref_dir1)",
        "staged_at(custom_plan.ref_dir2)",
        "staged_at(index)",
        "staged_at(custom_plan.marker_gene_table_ref1)",
        "staged_at(custom_plan.marker_gene_table_ref2)",
        "staged_at(custom_plan.pathway_map)",
        "staged_at(custom_plan.reaction_func_path)",
        "staged_at(custom_plan.regroup_map)",
    ):
        assert option in module

    assert "def add_asset = { value, flag ->" in main
    assert "custom_assets << source" in main
    assert "def has_mapping_overrides = [" in main
    assert "def add_override_asset = { value, flag ->" in main
    assert "override_assets << source" in main
    for staging_call in (
        "override_plan.pathway_map = add_override_asset",
        "override_plan.regroup_map = add_override_asset",
        "override_plan.reaction_func_path = add_override_asset",
        "override_plan.reaction_func_value = reaction_func.toString()",
    ):
        assert staging_call in main
    for option in (
        "params.picrust2_ref_dir1",
        "params.picrust2_ref_dir2",
        "params.picrust2_custom_trait_tables_ref1",
        "params.picrust2_custom_trait_tables_ref2",
        "params.picrust2_marker_gene_table_ref1",
        "params.picrust2_marker_gene_table_ref2",
        "params.picrust2_pathway_map",
        "params.picrust2_reaction_func",
        "params.picrust2_regroup_map",
    ):
        assert option in main
    for staging_call in (
        "custom_plan.ref_dir1 = add_asset(params.picrust2_ref_dir1",
        "values(params.picrust2_custom_trait_tables_ref1)",
        "add_asset(value, '--picrust2-custom-trait-tables-ref1')",
        "values(params.picrust2_custom_trait_tables_ref2)",
        "add_asset(value, '--picrust2-custom-trait-tables-ref2')",
        "custom_plan.marker_gene_table_ref1 = add_asset",
        "custom_plan.marker_gene_table_ref2 = add_asset",
        "custom_plan.pathway_map = add_optional_asset",
        "custom_plan.regroup_map = add_optional_asset",
        "custom_plan.reaction_func_path = add_asset",
    ):
        assert staging_call in main
    assert "MS_FUNCTIONAL_CUSTOM(" in main
    assert "MS_FUNCTIONAL_OVERRIDES(" in main
    assert "} else if (has_path_mapping_overrides) {" in main
    assert "} else if (has_mapping_overrides) {" in main
    assert "MS_FUNCTIONAL_SYMBOLIC_OVERRIDE(" in main
    assert "custom_assets_ch = Channel.value(custom_assets)" in main
    assert "custom_plan_ch = Channel.value(custom_plan)" in main
    assert "MS_FUNCTIONAL(MS_CLUSTER.out.table, MS_CLUSTER.out.rep_seqs)" in main
    assert (
        "def is_custom_database = params.picrust2_database?.toString()?.trim()?"
        ".equalsIgnoreCase('custom')"
    ) in main
    assert "if (has_custom_assets && !is_custom_database)" in main
    assert "if (is_custom_database)" in main


def test_nextflow_symbolic_reaction_override_does_not_bind_empty_path_input() -> None:
    module = (NEXTFLOW / "modules" / "ms_functional.nf").read_text(encoding="utf-8")
    main = (NEXTFLOW / "main.nf").read_text(encoding="utf-8")
    symbolic_process = module.split("process MS_FUNCTIONAL_SYMBOLIC_OVERRIDE", 1)[1].split(
        "// Custom references are a separate process", 1
    )[0]

    assert "path otu_table" in symbolic_process
    assert "path rep_seqs" in symbolic_process
    assert "override_assets" not in symbolic_process
    assert "path override_assets" not in symbolic_process
    assert "reaction_func_is_symbolic" in main
    assert "has_path_mapping_overrides" in main
    assert "MS_FUNCTIONAL_SYMBOLIC_OVERRIDE(" in main


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

    config = (NEXTFLOW / "nextflow.config").read_text(encoding="utf-8")

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
        # Images resolve from configurable registry/tag, not floating :ci tags.
        assert "${params.container_registry}" in text
        assert "${params.container_tag}" in text
        assert ":ci" not in text

    # Published GHCR registry is the default source; tag is pinnable.
    assert "ghcr.io/qwerty239qwe/microsuite" in config
    assert "container_registry" in config
    assert "container_tag" in config
    for image in [
        "fastqc",
        "multiqc",
        "qiime2-amplicon",
        "prjna321534-alpha",
        "microsuite-picrust2",
    ]:
        assert f"{image}:${{params.container_tag}}" in docker
    # Singularity pulls the same images over docker:// (no phantom .sif paths).
    assert "docker://${params.container_registry}" in singularity
    assert "containers/singularity" not in singularity
    assert ".sif'" not in singularity


def test_nextflow_docs_state_profiles_and_stub_status() -> None:
    docs = (ROOT / "docs" / "api-nextflow.md").read_text(encoding="utf-8")

    assert "-profile local" in docs
    assert "-profile docker" in docs
    assert "-profile singularity" in docs
    assert "module files contain runnable commands" in docs
    assert "Nextflow `-stub-run`" in docs


def test_nextflow_fastp_module_exists_and_declares_process() -> None:
    fastp = NEXTFLOW / "modules" / "fastp.nf"
    assert fastp.exists(), fastp
    text = fastp.read_text(encoding="utf-8")
    assert "process FASTP" in text
    assert "label 'fastp'" in text
    assert "params.fastp_cpus" in text
    assert "emit: trimmed" in text
    assert "emit: report" in text
    # PE and SE fastp invocations
    assert "--in1" in text and "--out1" in text
    assert "--in2" in text and "--out2" in text
    assert "params.fastp_args" in text
    assert "stub:" in text


def test_nextflow_fastp_params_and_container_labels() -> None:
    config = (NEXTFLOW / "nextflow.config").read_text(encoding="utf-8")
    for key in ("trim", "fastp_cpus", "fastp_args"):
        assert key in config, key
    docker = (NEXTFLOW / "profiles" / "docker.config").read_text(encoding="utf-8")
    singularity = (NEXTFLOW / "profiles" / "singularity.config").read_text(encoding="utf-8")
    assert "withLabel: fastp" in docker
    assert "withLabel: fastp" in singularity


def test_nextflow_main_wires_fastp_trim() -> None:
    main = (NEXTFLOW / "main.nf").read_text(encoding="utf-8")
    assert "include { FASTP } from './modules/fastp'" in main
    assert "FASTP(samples_ch)" in main
    assert "params.trim" in main
    assert "collectFile" in main
    assert "trimmed_manifest.tsv" in main
    assert ".mix(extra_qc)" in main
    # raw (non-breaking) path preserved in the else branch
    assert "dada2_manifest = manifest_ch" in main
    assert "dada2_reads  = reads_ch" in main or "dada2_reads = reads_ch" in main
    # FASTQC and DADA2 each still invoked exactly once, via the toggled inputs
    assert main.count("FASTQC(") == 1
    assert "QIIME2_DADA2(dada2_manifest, metadata_ch, dada2_reads)" in main
    # trimmed-manifest closure must normalize FASTP.out.trimmed's reads element
    # before indexing: single-end samples yield a scalar Path (not a List) from
    # the module's glob output, so unguarded reads.size()/reads[0] crashes.
    assert "reads instanceof List ? reads : [reads]" in main
