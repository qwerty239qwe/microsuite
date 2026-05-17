process QIIME2_TAXONOMY {
    tag 'qiime2_taxonomy'
    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path rep_seqs
    path classifier

    output:
    path 'taxonomy.qza', emit: taxonomy

    script:
    """
    echo "QIIME2 taxonomy placeholder for ${rep_seqs} and ${classifier}" > taxonomy.qza
    """
}
