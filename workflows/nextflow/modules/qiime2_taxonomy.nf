process QIIME2_TAXONOMY {
    tag 'qiime2_taxonomy'
    label 'qiime2'
    cpus { params.threads as int }
    publishDir "${params.outdir}/taxonomy", mode: 'copy'

    input:
    path rep_seqs
    path classifier

    output:
    path 'taxonomy.qza', emit: taxonomy

    script:
    """
    qiime feature-classifier classify-sklearn \
      --i-classifier ${classifier} \
      --i-reads ${rep_seqs} \
      --p-n-jobs ${task.cpus} \
      --o-classification taxonomy.qza
    """

    stub:
    """
    printf 'stub taxonomy\\n' > taxonomy.qza
    """
}
