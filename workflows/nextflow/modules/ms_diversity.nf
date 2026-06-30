process MS_DIVERSITY {
    tag 'ms_diversity'
    label 'microsuite_amplicon'
    publishDir "${params.outdir}/diversity", mode: 'copy'

    input:
    path table

    output:
    path 'alpha', emit: alpha_dir

    script:
    """
    mkdir -p alpha
    # breakaway and iNEXT are R-backed; the native metrics are pure Python.
    for metric in breakaway inext shannon observed_features chao1; do
      microsuite diversity alpha ${table} \
          --metric \${metric} \
          --output alpha/\${metric}.tsv \
        || echo "alpha \${metric} failed" >&2
    done
    """

    stub:
    """
    mkdir -p alpha
    printf 'sample_id\\tmethod\\testimate\\tstatus\\nS1\\tbreakaway\\t10\\tok\\n' > alpha/breakaway.tsv
    printf 'sample_id\\tshannon\\nS1\\t2.0\\n' > alpha/shannon.tsv
    """
}
