process REPORT {
    tag 'report'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path table
    path taxonomy
    path diversity_dir

    output:
    path 'report', emit: report_dir

    script:
    """
    mkdir -p report
    echo "Report placeholder for ${table}, ${taxonomy}, and ${diversity_dir}" > report/README.txt
    """
}
