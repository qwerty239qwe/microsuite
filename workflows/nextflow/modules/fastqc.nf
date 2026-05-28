process FASTQC {
    tag "${sample_id}"
    label 'fastqc'
    cpus 2
    publishDir "${params.outdir}/qc/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path 'fastqc', emit: qc_dir

    script:
    def read_args = reads instanceof List ? reads.collect { it.toString() }.join(' ') : reads.toString()
    """
    mkdir -p fastqc
    fastqc --outdir fastqc --threads ${task.cpus} ${read_args}
    """

    stub:
    def read_args = reads instanceof List ? reads.collect { it.toString() }.join(' ') : reads.toString()
    """
    mkdir -p fastqc
    for read in ${read_args}; do
      name=\$(basename "\${read}")
      base="\${name%%.*}"
      mkdir -p "fastqc/\${base}_fastqc"
      printf '##FastQC\\tstub\\n' > "fastqc/\${base}_fastqc/fastqc_data.txt"
      printf '<html><body>stub FastQC for %s</body></html>\\n' "\${name}" > "fastqc/\${base}_fastqc.html"
      printf 'stub zip for %s\\n' "\${name}" > "fastqc/\${base}_fastqc.zip"
    done
    """
}
