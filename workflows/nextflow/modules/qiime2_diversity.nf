process QIIME2_DIVERSITY {
    tag 'qiime2_diversity'
    publishDir "${params.outdir}/diversity", mode: 'copy'

    input:
    path table
    path rooted_tree

    output:
    path 'diversity', emit: diversity_dir

    script:
    """
    mkdir -p diversity
    echo "QIIME2 diversity placeholder for ${table} and ${rooted_tree}" > diversity/README.txt
    """
}
