process QIIME2_DADA2 {
    tag 'qiime2_dada2'
    publishDir "${params.outdir}/denoise", mode: 'copy'

    input:
    path manifest
    path metadata

    output:
    path 'table.qza', emit: table
    path 'rep-seqs.qza', emit: rep_seqs

    script:
    """
    echo "QIIME2 DADA2 placeholder for ${manifest} and ${metadata}" > table.qza
    echo "QIIME2 DADA2 representative sequences placeholder" > rep-seqs.qza
    """
}
