process MS_FUNCTIONAL {
    tag 'ms_functional'
    label 'microsuite_picrust2'
    cpus { params.threads as int }
    publishDir "${params.outdir}/functional", mode: 'copy'

    input:
    path otu_table
    path rep_seqs

    output:
    path 'picrust2', emit: functional_dir

    script:
    """
    # PICRUSt2 predicts gene-family/pathway abundance from 16S OTUs + counts.
    microsuite functional_profile \
        --backend picrust2 \
        --table ${otu_table} \
        --rep-seqs ${rep_seqs} \
        --output-dir picrust2 \
        --threads ${task.cpus}
    """

    stub:
    """
    mkdir -p picrust2
    printf 'stub picrust2 output\\n' > picrust2/pathways.tsv
    """
}
