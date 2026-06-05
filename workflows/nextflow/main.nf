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
params.threads = params.threads ?: 2
params.trim_left = params.trim_left ?: 0
params.trunc_len = params.trunc_len ?: 0
params.trim_left_f = params.trim_left_f ?: 0
params.trunc_len_f = params.trunc_len_f ?: 0
params.trim_left_r = params.trim_left_r ?: 0
params.trunc_len_r = params.trunc_len_r ?: 0
params.sampling_depth = params.sampling_depth ?: 1000

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

    FASTQC(samples_ch)
    MULTIQC(FASTQC.out.qc_dir.collect())
    QIIME2_DADA2(manifest_ch, metadata_ch, reads_ch)
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
