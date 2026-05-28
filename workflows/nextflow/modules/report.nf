process REPORT {
    tag 'report'
    label 'microsuite'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path table
    path taxonomy
    path diversity_dir
    path multiqc_dir

    output:
    path 'report.html', emit: html
    path 'run.json', emit: run_json

    script:
    """
    cat > report.html <<'EOF'
    <!doctype html>
    <html lang="en">
    <head><meta charset="utf-8"><title>microsuite amplicon_qiime2 report</title></head>
    <body>
      <h1>microsuite amplicon_qiime2</h1>
      <ul>
        <li>Feature table: ${table}</li>
        <li>Taxonomy: ${taxonomy}</li>
        <li>Diversity directory: ${diversity_dir}</li>
        <li>MultiQC directory: ${multiqc_dir}</li>
      </ul>
    </body>
    </html>
    EOF

    cat > run.json <<'EOF'
    {
      "workflow": "amplicon_qiime2",
      "table": "${table}",
      "taxonomy": "${taxonomy}",
      "diversity": "${diversity_dir}",
      "multiqc": "${multiqc_dir}"
    }
    EOF
    """

    stub:
    """
    printf '<html><body>stub report</body></html>\\n' > report.html
    printf '{"workflow":"amplicon_qiime2","stub":true}\\n' > run.json
    """
}
