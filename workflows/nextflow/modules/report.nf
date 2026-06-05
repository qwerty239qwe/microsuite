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
    {
      printf '%s\\n' '<!doctype html>'
      printf '%s\\n' '<html lang="en">'
      printf '%s\\n' '<head><meta charset="utf-8"><title>microsuite amplicon_qiime2 report</title></head>'
      printf '%s\\n' '<body>'
      printf '%s\\n' '<h1>microsuite amplicon_qiime2</h1>'
      printf '%s\\n' '<ul>'
      printf '%s\\n' '<li>Feature table: ${table}</li>'
      printf '%s\\n' '<li>Taxonomy: ${taxonomy}</li>'
      printf '%s\\n' '<li>Diversity directory: ${diversity_dir}</li>'
      printf '%s\\n' '<li>MultiQC directory: ${multiqc_dir}</li>'
      printf '%s\\n' '</ul>'
      printf '%s\\n' '</body>'
      printf '%s\\n' '</html>'
    } > report.html

    {
      printf '%s\\n' '{'
      printf '%s\\n' '  "workflow": "amplicon_qiime2",'
      printf '%s\\n' '  "table": "${table}",'
      printf '%s\\n' '  "taxonomy": "${taxonomy}",'
      printf '%s\\n' '  "diversity": "${diversity_dir}",'
      printf '%s\\n' '  "multiqc": "${multiqc_dir}"'
      printf '%s\\n' '}'
    } > run.json
    """

    stub:
    """
    printf '<html><body>stub report</body></html>\\n' > report.html
    printf '{"workflow":"amplicon_qiime2","stub":true}\\n' > run.json
    """
}
