nextflow.enable.dsl = 2

include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'
include { QIIME2_DADA2 } from './modules/qiime2_dada2'
include { QIIME2_TAXONOMY } from './modules/qiime2_taxonomy'
include { QIIME2_PHYLOGENY } from './modules/qiime2_phylogeny'
include { QIIME2_DIVERSITY } from './modules/qiime2_diversity'
include { REPORT } from './modules/report'

params.workflow = params.workflow ?: 'amplicon_qiime2'
params.manifest = params.manifest ?: null
params.metadata = params.metadata ?: null
params.classifier = params.classifier ?: null
params.outdir = params.outdir ?: 'results'

workflow {
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

    manifest_ch = Channel.value(file(params.manifest))
    metadata_ch = Channel.value(file(params.metadata))
    classifier_ch = Channel.value(file(params.classifier))

    FASTQC(manifest_ch)
    MULTIQC(FASTQC.out.qc_dir)
    QIIME2_DADA2(manifest_ch, metadata_ch)
    QIIME2_TAXONOMY(QIIME2_DADA2.out.rep_seqs, classifier_ch)
    QIIME2_PHYLOGENY(QIIME2_DADA2.out.rep_seqs)
    QIIME2_DIVERSITY(QIIME2_DADA2.out.table, QIIME2_PHYLOGENY.out.rooted_tree)
    REPORT(QIIME2_DADA2.out.table, QIIME2_TAXONOMY.out.taxonomy, QIIME2_DIVERSITY.out.diversity_dir)
}
