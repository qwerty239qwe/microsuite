process QIIME2_PHYLOGENY {
    tag 'qiime2_phylogeny'
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path rep_seqs

    output:
    path 'rooted-tree.qza', emit: rooted_tree

    script:
    """
    echo "QIIME2 phylogeny placeholder for ${rep_seqs}" > rooted-tree.qza
    """
}
