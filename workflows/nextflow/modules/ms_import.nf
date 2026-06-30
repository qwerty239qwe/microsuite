process MS_IMPORT {
    tag 'ms_import'
    label 'microsuite_amplicon'
    publishDir "${params.outdir}/import", mode: 'copy'

    input:
    path otu_table
    path metadata

    output:
    path 'table.h5ad', emit: table

    script:
    """
    microsuite import tsv ${otu_table} \
        --metadata ${metadata} \
        --output table.h5ad
    """

    stub:
    """
    printf 'stub h5ad\\n' > table.h5ad
    """
}
