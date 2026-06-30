process MS_CLUSTER {
    tag 'ms_cluster'
    label 'microsuite_amplicon'
    cpus { params.threads as int }
    publishDir "${params.outdir}/cluster", mode: 'copy'

    input:
    path reads_fasta

    output:
    path 'otu_table.tsv', emit: table
    path 'otus.fasta', emit: rep_seqs

    script:
    """
    # Dereplicate with singleton removal first (microsuite has no derep CLI),
    # then drive OTU picking + per-sample counting through the microsuite CLI.
    vsearch --derep_fulllength ${reads_fasta} \
        --output uniques.fasta \
        --sizeout \
        --minuniquesize 2 \
        --threads ${task.cpus}

    microsuite cluster \
        --backend vsearch \
        --rep-seqs uniques.fasta \
        --reads ${reads_fasta} \
        --output-table otu_table.tsv \
        --output-rep-seqs otus.fasta \
        --identity ${params.otu_identity} \
        --sample-delimiter '_' \
        --sample-field 0
    """

    stub:
    """
    printf 'feature-id\\tS1\\tS2\\nOTU_1\\t5\\t3\\n' > otu_table.tsv
    printf '>OTU_1\\nACGT\\n' > otus.fasta
    """
}
