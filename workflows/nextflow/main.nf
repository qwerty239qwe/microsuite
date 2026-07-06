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
    MS_FUNCTIONAL(MS_CLUSTER.out.table, MS_CLUSTER.out.rep_seqs)
    MS_REPORT(
        MS_CLUSTER.out.table,
        MS_DIVERSITY.out.alpha_dir,
        MS_FUNCTIONAL.out.functional_dir
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
