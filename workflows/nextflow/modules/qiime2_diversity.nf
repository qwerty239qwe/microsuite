process QIIME2_DIVERSITY {
    tag 'qiime2_diversity'
    label 'qiime2'
    publishDir "${params.outdir}/diversity", mode: 'copy'

    input:
    path table
    path rooted_tree
    path metadata

    output:
    path 'diversity', emit: diversity_dir

    script:
    """
    qiime diversity core-metrics-phylogenetic \
      --i-table ${table} \
      --i-phylogeny ${rooted_tree} \
      --m-metadata-file ${metadata} \
      --p-sampling-depth ${params.sampling_depth} \
      --output-dir diversity
    """

    stub:
    """
    mkdir -p diversity
    printf 'stub core metrics\\n' > diversity/README.txt
    printf 'stub emperor\\n' > diversity/bray_curtis_emperor.qzv
    """
}
