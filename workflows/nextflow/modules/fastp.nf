process FASTP {
    tag "${sample_id}"
    label 'fastp'
    cpus { params.fastp_cpus as int }
    publishDir "${params.outdir}/trim/fastp", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}*.trim.fastq.gz"), emit: trimmed
    path "${sample_id}.fastp.json", emit: report
    path "${sample_id}.fastp.html"

    script:
    def paired = reads instanceof List && reads.size() > 1
    if (paired)
        """
        fastp --in1 ${reads[0]} --in2 ${reads[1]} \
              --out1 ${sample_id}_1.trim.fastq.gz --out2 ${sample_id}_2.trim.fastq.gz \
              --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html \
              --thread ${task.cpus} ${params.fastp_args}
        """
    else
        """
        fastp --in1 ${reads instanceof List ? reads[0] : reads} \
              --out1 ${sample_id}.trim.fastq.gz \
              --json ${sample_id}.fastp.json --html ${sample_id}.fastp.html \
              --thread ${task.cpus} ${params.fastp_args}
        """

    stub:
    def paired = reads instanceof List && reads.size() > 1
    if (paired)
        """
        printf '' | gzip > ${sample_id}_1.trim.fastq.gz
        printf '' | gzip > ${sample_id}_2.trim.fastq.gz
        printf '{"summary":{}}\\n' > ${sample_id}.fastp.json
        printf '<html><body>stub fastp</body></html>\\n' > ${sample_id}.fastp.html
        """
    else
        """
        printf '' | gzip > ${sample_id}.trim.fastq.gz
        printf '{"summary":{}}\\n' > ${sample_id}.fastp.json
        printf '<html><body>stub fastp</body></html>\\n' > ${sample_id}.fastp.html
        """
}
