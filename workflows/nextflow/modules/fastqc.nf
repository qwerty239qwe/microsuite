process FASTQC {
    tag 'fastqc'
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    path manifest

    output:
    path 'fastqc', emit: qc_dir

    script:
    """
    mkdir -p fastqc
    echo "FASTQC placeholder for ${manifest}" > fastqc/README.txt
    """
}
