process QIIME2_DADA2 {
    tag 'qiime2_dada2'
    label 'qiime2'
    cpus { params.threads as int }
    publishDir "${params.outdir}/denoise", mode: 'copy'

    input:
    path manifest
    path metadata
    path reads

    output:
    path 'table.qza', emit: table
    path 'rep-seqs.qza', emit: rep_seqs
    path 'stats.qza', emit: stats

    script:
    """
    awk -F '\\t' '
      BEGIN { OFS="\\t"; print "sample-id", "absolute-filepath", "direction" }
      NR == 1 { for (i = 1; i <= NF; i++) h[\$i] = i; next }
      {
        sample = \$h["sample_id"]
        r1 = \$h["read1"]
        gsub(/^.*[\\\\/]/, "", r1)
        print sample, ENVIRON["PWD"] "/" r1, "forward"
        if (("read2" in h) && \$h["read2"] != "") {
          r2 = \$h["read2"]
          gsub(/^.*[\\\\/]/, "", r2)
          print sample, ENVIRON["PWD"] "/" r2, "reverse"
        }
      }
    ' ${manifest} > qiime-manifest.tsv

    paired=\$(awk -F '\\t' '
      NR == 1 { for (i = 1; i <= NF; i++) h[\$i] = i; next }
      ("read2" in h) && \$h["read2"] != "" { found = 1 }
      END { print found ? "true" : "false" }
    ' ${manifest})

    if [ "\${paired}" = "true" ]; then
      qiime tools import \
        --type 'SampleData[PairedEndSequencesWithQuality]' \
        --input-path qiime-manifest.tsv \
        --output-path demux.qza \
        --input-format PairedEndFastqManifestPhred33V2
      qiime dada2 denoise-paired \
        --i-demultiplexed-seqs demux.qza \
        --p-trim-left-f ${params.trim_left_f} \
        --p-trunc-len-f ${params.trunc_len_f} \
        --p-trim-left-r ${params.trim_left_r} \
        --p-trunc-len-r ${params.trunc_len_r} \
        --p-n-threads ${task.cpus} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats stats.qza
    else
      qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path qiime-manifest.tsv \
        --output-path demux.qza \
        --input-format SingleEndFastqManifestPhred33V2
      qiime dada2 denoise-single \
        --i-demultiplexed-seqs demux.qza \
        --p-trim-left ${params.trim_left} \
        --p-trunc-len ${params.trunc_len} \
        --p-n-threads ${task.cpus} \
        --o-table table.qza \
        --o-representative-sequences rep-seqs.qza \
        --o-denoising-stats stats.qza
    fi
    """

    stub:
    """
    printf 'stub table\\n' > table.qza
    printf 'stub representative sequences\\n' > rep-seqs.qza
    printf 'stub denoising stats\\n' > stats.qza
    """
}
