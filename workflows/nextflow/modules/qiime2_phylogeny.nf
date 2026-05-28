process QIIME2_PHYLOGENY {
    tag 'qiime2_phylogeny'
    label 'qiime2'
    publishDir "${params.outdir}/phylogeny", mode: 'copy'

    input:
    path rep_seqs

    output:
    path 'aligned-rep-seqs.qza', emit: aligned_rep_seqs
    path 'masked-aligned-rep-seqs.qza', emit: masked_aligned_rep_seqs
    path 'unrooted-tree.qza', emit: unrooted_tree
    path 'rooted-tree.qza', emit: rooted_tree

    script:
    """
    qiime phylogeny align-to-tree-mafft-fasttree \
      --i-sequences ${rep_seqs} \
      --o-alignment aligned-rep-seqs.qza \
      --o-masked-alignment masked-aligned-rep-seqs.qza \
      --o-tree unrooted-tree.qza \
      --o-rooted-tree rooted-tree.qza
    """

    stub:
    """
    printf 'stub alignment\\n' > aligned-rep-seqs.qza
    printf 'stub masked alignment\\n' > masked-aligned-rep-seqs.qza
    printf 'stub unrooted tree\\n' > unrooted-tree.qza
    printf 'stub rooted tree\\n' > rooted-tree.qza
    """
}
