nextflow.enable.dsl = 2

include { FASTQC } from './modules/fastqc'
include { FASTP } from './modules/fastp'
include { MULTIQC } from './modules/multiqc'
include { QIIME2_DADA2 } from './modules/qiime2_dada2'
include { QIIME2_TAXONOMY } from './modules/qiime2_taxonomy'
include { QIIME2_PHYLOGENY } from './modules/qiime2_phylogeny'
include { QIIME2_DIVERSITY } from './modules/qiime2_diversity'
include { REPORT } from './modules/report'
include { MS_CLUSTER } from './modules/ms_cluster'
include { MS_IMPORT } from './modules/ms_import'
include { MS_DIVERSITY } from './modules/ms_diversity'
include { MS_FUNCTIONAL } from './modules/ms_functional'
include { MS_FUNCTIONAL_OVERRIDES } from './modules/ms_functional'
include { MS_FUNCTIONAL_SYMBOLIC_OVERRIDE } from './modules/ms_functional'
include { MS_FUNCTIONAL_CUSTOM } from './modules/ms_functional'
include { MS_REPORT } from './modules/ms_report'

params.workflow = params.workflow ?: 'amplicon_qiime2'
params.manifest = params.manifest ?: null
params.metadata = params.metadata ?: null
params.classifier = params.classifier ?: null
params.outdir = params.outdir ?: 'results'
params.threads = params.threads ?: 2
params.trim_left = params.trim_left ?: 0
params.trunc_len = params.trunc_len ?: 0
params.trim_left_f = params.trim_left_f ?: 0
params.trunc_len_f = params.trunc_len_f ?: 0
params.trim_left_r = params.trim_left_r ?: 0
params.trunc_len_r = params.trunc_len_r ?: 0
params.sampling_depth = params.sampling_depth ?: 1000
// amplicon_microsuite inputs: a combined FASTA of per-sample, sample-labelled
// reads (e.g. ">SRR..._1;sample=SRR...;") plus sample metadata.
params.reads_fasta = params.reads_fasta ?: null
params.otu_identity = params.otu_identity ?: 0.97

def resolveManifestPath(manifest_path, raw_path) {
    if (raw_path == null) {
        return null
    }
    def value = raw_path.toString().trim()
    if (!value) {
        return null
    }
    def candidate = file(value)
    return candidate.isAbsolute() ? candidate : file("${manifest_path.parent}/${value}")
}

workflow amplicon_microsuite {
    // Drives the microsuite CLI end to end on a combined, sample-labelled read
    // FASTA: cluster -> import -> diversity + functional profiling -> report.
    if (!params.reads_fasta) {
        error "Missing required parameter: --reads_fasta"
    }
    if (!params.metadata) {
        error "Missing required parameter: --metadata"
    }

    reads_ch = Channel.value(file(params.reads_fasta))
    metadata_ch = Channel.value(file(params.metadata))

    MS_CLUSTER(reads_ch)
    MS_IMPORT(MS_CLUSTER.out.table, metadata_ch)
    MS_DIVERSITY(MS_IMPORT.out.table)
    def functional_dir
    def is_custom_database = params.picrust2_database?.toString()?.trim()?.equalsIgnoreCase('custom')
    def has_custom_asset_param = { value ->
        value instanceof Collection
            ? !value.isEmpty()
            : value != null && value.toString().trim()
    }
    def has_reaction_func = has_custom_asset_param(params.picrust2_reaction_func)
    def reaction_func_is_path = has_reaction_func && file(params.picrust2_reaction_func.toString()).exists()
    def reaction_func_is_symbolic = has_reaction_func &&
        (params.picrust2_reaction_func.toString().trim() ==~ /[A-Za-z][A-Za-z0-9_.-]*/)
    if (has_reaction_func && !reaction_func_is_path && !reaction_func_is_symbolic) {
        error '--picrust2-reaction-func must be an existing path or a symbolic value such as EC'
    }
    def has_custom_assets = [
        params.picrust2_ref_dir1,
        params.picrust2_ref_dir2,
        params.picrust2_custom_trait_tables_ref1,
        params.picrust2_custom_trait_tables_ref2,
        params.picrust2_marker_gene_table_ref1,
        params.picrust2_marker_gene_table_ref2,
    ].any { value -> has_custom_asset_param(value) }
    def has_mapping_overrides = [
        params.picrust2_pathway_map,
        params.picrust2_reaction_func,
        params.picrust2_regroup_map,
    ].any { value -> has_custom_asset_param(value) }
    def has_path_mapping_overrides = [
        params.picrust2_pathway_map,
        params.picrust2_regroup_map,
    ].any { value -> has_custom_asset_param(value) } || reaction_func_is_path
    if (has_custom_assets && !is_custom_database) {
        error 'PICRUSt2 custom reference options require --picrust2_database custom'
    }
    if (is_custom_database) {
        // Build one path collection containing every custom asset. Nextflow
        // stages this collection into the task directory; custom_plan stores
        // only indexes into that staged collection, never host path strings.
        def custom_assets = []
        def custom_plan = [
            trait_tables_ref1: [],
            trait_tables_ref2: [],
            reaction_func_path: null,
            reaction_func_value: null,
        ]
        def values = { value ->
            if (value == null) {
                return []
            }
            if (value instanceof Collection) {
                return value.findAll { it != null && it.toString().trim() }
            }
            return value.toString().split(',').collect { it.trim() }.findAll { it }
        }
        def add_asset = { value, flag ->
            if (value == null || !value.toString().trim()) {
                error "${flag} is required for --picrust2-database custom"
            }
            def source = file(value.toString())
            if (!source.exists()) {
                error "${flag} does not exist: ${value}"
            }
            custom_assets << source
            custom_assets.size() - 1
        }
        def add_optional_asset = { value, flag ->
            value == null || !value.toString().trim() ? null : add_asset(value, flag)
        }

        custom_plan.ref_dir1 = add_asset(params.picrust2_ref_dir1, '--picrust2-ref-dir1')
        values(params.picrust2_custom_trait_tables_ref1).each { value ->
            custom_plan.trait_tables_ref1 << add_asset(value, '--picrust2-custom-trait-tables-ref1')
        }
        if (custom_plan.trait_tables_ref1.isEmpty()) {
            error '--picrust2-custom-trait-tables-ref1 is required for --picrust2-database custom'
        }
        custom_plan.marker_gene_table_ref1 = add_asset(
            params.picrust2_marker_gene_table_ref1,
            '--picrust2-marker-gene-table-ref1'
        )

        def ref2_values = values(params.picrust2_custom_trait_tables_ref2)
        def has_ref2 = params.picrust2_ref_dir2 != null || !ref2_values.isEmpty() ||
            (params.picrust2_marker_gene_table_ref2 != null && params.picrust2_marker_gene_table_ref2.toString().trim())
        if (has_ref2) {
            custom_plan.ref_dir2 = add_asset(params.picrust2_ref_dir2, '--picrust2-ref-dir2')
            ref2_values.each { value ->
                custom_plan.trait_tables_ref2 << add_asset(value, '--picrust2-custom-trait-tables-ref2')
            }
            if (custom_plan.trait_tables_ref2.isEmpty()) {
                error '--picrust2-custom-trait-tables-ref2 is required when reference 2 is configured'
            }
            custom_plan.marker_gene_table_ref2 = add_asset(
                params.picrust2_marker_gene_table_ref2,
                '--picrust2-marker-gene-table-ref2'
            )
        } else {
            custom_plan.ref_dir2 = null
            custom_plan.marker_gene_table_ref2 = null
        }

        custom_plan.pathway_map = add_optional_asset(params.picrust2_pathway_map, '--picrust2-pathway-map')
        custom_plan.regroup_map = add_optional_asset(params.picrust2_regroup_map, '--picrust2-regroup-map')
        def reaction_func = params.picrust2_reaction_func
        if (reaction_func != null && reaction_func.toString().trim()) {
            def reaction_candidate = file(reaction_func.toString())
            if (reaction_candidate.exists()) {
                custom_plan.reaction_func_path = add_asset(
                    reaction_func,
                    '--picrust2-reaction-func'
                )
            } else {
                // Upstream also accepts a symbolic reaction/function name.
                custom_plan.reaction_func_value = reaction_func.toString()
            }
        }

        def custom_assets_ch = Channel.value(custom_assets)
        def custom_plan_ch = Channel.value(custom_plan)
        MS_FUNCTIONAL_CUSTOM(
            MS_CLUSTER.out.table,
            MS_CLUSTER.out.rep_seqs,
            custom_assets_ch,
            custom_plan_ch
        )
        functional_dir = MS_FUNCTIONAL_CUSTOM.out.functional_dir
    } else if (has_path_mapping_overrides) {
        // SC/oldIMG accept mapping overrides without custom references. Stage
        // every path-valued override and keep symbolic reaction functions as
        // values in the non-path plan.
        def override_assets = []
        def override_plan = [
            pathway_map: null,
            reaction_func_path: null,
            reaction_func_value: null,
            regroup_map: null,
        ]
        def add_override_asset = { value, flag ->
            def source = file(value.toString())
            if (!source.exists()) {
                error "${flag} does not exist: ${value}"
            }
            override_assets << source
            override_assets.size() - 1
        }
        if (params.picrust2_pathway_map != null && params.picrust2_pathway_map.toString().trim()) {
            override_plan.pathway_map = add_override_asset(
                params.picrust2_pathway_map,
                '--picrust2-pathway-map'
            )
        }
        if (params.picrust2_regroup_map != null && params.picrust2_regroup_map.toString().trim()) {
            override_plan.regroup_map = add_override_asset(
                params.picrust2_regroup_map,
                '--picrust2-regroup-map'
            )
        }
        def reaction_func = params.picrust2_reaction_func
        if (reaction_func != null && reaction_func.toString().trim()) {
            def reaction_candidate = file(reaction_func.toString())
            if (reaction_candidate.exists()) {
                override_plan.reaction_func_path = add_override_asset(
                    reaction_func,
                    '--picrust2-reaction-func'
                )
            } else {
                override_plan.reaction_func_value = reaction_func.toString()
            }
        }

        def override_assets_ch = Channel.value(override_assets)
        def override_plan_ch = Channel.value(override_plan)
        MS_FUNCTIONAL_OVERRIDES(
            MS_CLUSTER.out.table,
            MS_CLUSTER.out.rep_seqs,
            override_assets_ch,
            override_plan_ch
        )
        functional_dir = MS_FUNCTIONAL_OVERRIDES.out.functional_dir
    } else if (has_mapping_overrides) {
        // A symbolic reaction function has no asset to stage. Use the
        // process with only the table and representative-sequence path inputs.
        MS_FUNCTIONAL_SYMBOLIC_OVERRIDE(
            MS_CLUSTER.out.table,
            MS_CLUSTER.out.rep_seqs
        )
        functional_dir = MS_FUNCTIONAL_SYMBOLIC_OVERRIDE.out.functional_dir
    } else {
        // Preserve the original two-input call for SC and oldIMG.
        MS_FUNCTIONAL(MS_CLUSTER.out.table, MS_CLUSTER.out.rep_seqs)
        functional_dir = MS_FUNCTIONAL.out.functional_dir
    }
    MS_REPORT(
        MS_CLUSTER.out.table,
        MS_DIVERSITY.out.alpha_dir,
        functional_dir
    )
}

workflow {
    if (params.workflow == 'amplicon_microsuite') {
        amplicon_microsuite()
        return
    }
    if (params.workflow != 'amplicon_qiime2') {
        error "Unsupported workflow: ${params.workflow}"
    }
    if (!params.manifest) {
        error "Missing required parameter: --manifest"
    }
    if (!params.metadata) {
        error "Missing required parameter: --metadata"
    }
    if (!params.classifier) {
        error "Missing required parameter: --classifier"
    }

    manifest_path = file(params.manifest)
    manifest_ch = Channel.value(manifest_path)
    metadata_ch = Channel.value(file(params.metadata))
    classifier_ch = Channel.value(file(params.classifier))

    samples_ch = Channel
        .fromPath(params.manifest)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def sample_id = row.sample_id?.toString()
            if (!sample_id) {
                error "Manifest rows must include sample_id"
            }
            def read1 = resolveManifestPath(manifest_path, row.read1)
            if (read1 == null) {
                error "Manifest row ${sample_id} is missing read1"
            }
            def reads = [read1]
            def read2 = resolveManifestPath(manifest_path, row.read2)
            if (read2 != null) {
                reads << read2
            }
            tuple(sample_id, reads)
        }
    reads_ch = samples_ch.map { sample_id, reads -> reads }.flatten().collect()

    if (params.trim) {
        FASTP(samples_ch)
        trimmed_ch = FASTP.out.trimmed
        dada2_manifest = trimmed_ch
            .map { sid, reads ->
                def rs = reads instanceof List ? reads : [reads]
                def r2 = rs.size() > 1 ? rs[1].name : ''
                "${sid}\t${rs[0].name}\t${r2}\n"
            }
            .collectFile(name: 'trimmed_manifest.tsv', seed: 'sample_id\tread1\tread2\n', sort: true)
        fastqc_input = trimmed_ch
        dada2_reads  = trimmed_ch.map { sid, reads -> reads }.flatten().collect()
        extra_qc     = FASTP.out.report
    } else {
        fastqc_input = samples_ch
        dada2_manifest = manifest_ch
        dada2_reads  = reads_ch
        extra_qc     = Channel.empty()
    }

    FASTQC(fastqc_input)
    MULTIQC(FASTQC.out.qc_dir.mix(extra_qc).collect())
    QIIME2_DADA2(dada2_manifest, metadata_ch, dada2_reads)
    QIIME2_TAXONOMY(QIIME2_DADA2.out.rep_seqs, classifier_ch)
    QIIME2_PHYLOGENY(QIIME2_DADA2.out.rep_seqs)
    QIIME2_DIVERSITY(QIIME2_DADA2.out.table, QIIME2_PHYLOGENY.out.rooted_tree, metadata_ch)
    REPORT(
        QIIME2_DADA2.out.table,
        QIIME2_TAXONOMY.out.taxonomy,
        QIIME2_DIVERSITY.out.diversity_dir,
        MULTIQC.out.report_dir
    )
}
