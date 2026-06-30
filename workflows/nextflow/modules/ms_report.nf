process MS_REPORT {
    tag 'ms_report'
    label 'microsuite_amplicon'
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
    path otu_table
    path alpha_dir
    path functional_dir

    output:
    path 'report.html', emit: html
    path 'run.json', emit: run_json

    script:
    """
    {
      printf '%s\\n' '<!doctype html>'
      printf '%s\\n' '<html lang="en">'
      printf '%s\\n' '<head><meta charset="utf-8"><title>microsuite amplicon_microsuite report</title></head>'
      printf '%s\\n' '<body>'
      printf '%s\\n' '<h1>microsuite amplicon_microsuite</h1>'
      printf '%s\\n' '<ul>'
      printf '%s\\n' '<li>OTU table: ${otu_table}</li>'
      printf '%s\\n' '<li>Alpha diversity: ${alpha_dir}</li>'
      printf '%s\\n' '<li>Functional profile: ${functional_dir}</li>'
      printf '%s\\n' '</ul>'
      printf '%s\\n' '</body>'
      printf '%s\\n' '</html>'
    } > report.html

    {
      printf '%s\\n' '{'
      printf '%s\\n' '  "workflow": "amplicon_microsuite",'
      printf '%s\\n' '  "otu_table": "${otu_table}",'
      printf '%s\\n' '  "alpha": "${alpha_dir}",'
      printf '%s\\n' '  "functional": "${functional_dir}"'
      printf '%s\\n' '}'
    } > run.json
    """

    stub:
    """
    printf '<html><body>amplicon_microsuite</body></html>\\n' > report.html
    printf '{"workflow": "amplicon_microsuite"}\\n' > run.json
    """
}
