process MULTIQC {
    tag 'multiqc'
    label 'multiqc'
    publishDir "${params.outdir}/qc/multiqc", mode: 'copy'

    input:
    path qc_dirs

    output:
    path 'multiqc', emit: report_dir

    script:
    """
    mkdir -p multiqc
    multiqc ${qc_dirs.join(' ')} --outdir multiqc --force
    """

    stub:
    """
    mkdir -p multiqc
    printf '<html><body>stub MultiQC</body></html>\\n' > multiqc/multiqc_report.html
    """
}
