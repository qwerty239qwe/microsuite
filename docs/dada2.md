# DADA2 denoising

microsuite targets DADA2 `1.40.0` on Bioconductor `3.23` for direct R workflows
and q2-dada2 `2026.4.0` for QIIME 2 artifact workflows. The wrapper preserves
the legacy default: single-end DADA2 unless `paired=True`, `--paired`, or
`mode="paired"` is supplied.

## Backend choice

Use `qiime2-dada2` when the input is already a QIIME 2 demultiplexed reads
artifact and the outputs should remain QIIME artifacts. This backend exposes
q2-dada2 `denoise-single`, `denoise-paired`, `denoise-ccs`, and
`denoise-pyro`; PacBio CCS data should use `mode="ccs"`.

Use `dada2-r` when the input is a FASTQ directory and the desired outputs are
microsuite-native files: an ASV table TSV, representative sequence FASTA, and
denoising stats TSV.

## Primers, ITS, and multi-run studies

For ITS workflows, remove primers first with an existing trimming backend such
as Cutadapt or Trim Galore, then run DADA2 on the trimmed reads. Do not rely on
fixed-length truncation alone for variable-length ITS amplicons.

Large multi-run studies are not fully automated yet. The recommended pattern is
to denoise each sequencing run separately with identical parameters where
possible, merge the sequence tables, then remove chimeras and assign taxonomy
globally.

## Diagnostics

For q2-dada2 runs, `--output-base-transition-stats` writes the DADA2 base
transition stats artifact. Supplying both `--output-base-transition-stats` and
`--output-base-transition-plot` also runs `qiime dada2 plot-base-transitions`
and writes a `.qzv` visualization.
