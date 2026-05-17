process MULTIQC {
    tag 'multiqc'
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path qc_dir

    output:
    path 'multiqc', emit: report_dir

    script:
    """
    mkdir -p multiqc
    echo "MULTIQC placeholder for ${qc_dir}" > multiqc/README.txt
    """
}
