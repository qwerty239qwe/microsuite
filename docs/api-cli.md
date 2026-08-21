# CLI API

The CLI API is for local, one-step microbiome tasks.

Use it when you want to run a single operation, debug a workflow step, or teach a
small example:

```bash
microsuite qc --backend fastqc --input sample_R1.fastq.gz --output-dir qc/fastqc
microsuite qc --backend multiqc --input-dir qc/fastqc --output-dir qc/multiqc
microsuite qc --backend qiime2-demux --demux demux.qza -o qc/demux.qzv
microsuite trim --backend fastp --read1 sample_R1.fastq.gz --output1 trimmed_R1.fastq.gz --html qc/fastp.html --json-report qc/fastp.json
microsuite primer-check --input sample_R1.fastq.gz --front '^TCAGNNNNNNNNNNGGATTAGATACCCTGGTAGT' -o qc/primer-check.json
microsuite denoise --backend qiime2-dada2 --demux demux.qza --output-table table.qza --output-rep-seqs rep-seqs.qza --output-stats stats.qza --trunc-len 150
microsuite cluster --backend vsearch --table table.qza --rep-seqs rep-seqs.qza --output-table clustered-table.qza --output-rep-seqs clustered-rep-seqs.qza
microsuite normalize --backend native --method clr --table table.h5ad -o clr.h5ad
microsuite abundance --backend native --table table.h5ad --level genus -o abundance.tsv
microsuite tax_classify --backend qiime2 --rep-seqs rep-seqs.qza --classifier classifier.qza -o taxonomy.qza
microsuite diff_abundance --backend ancombc --table table.h5ad --group treatment -o diff.tsv
microsuite diversity beta-significance beta.tsv --metadata metadata.tsv --backend vegan --formula "site + phase" --strata subject --method adonis2 --runtime docker -o adonis2.tsv
microsuite diversity beta-significance beta.tsv --metadata metadata.tsv --backend vegan --formula "batch + group / time_point" --strata "batch:subject" --permutations 999 --method adonis2 --runtime docker -o repeated-measures.tsv
microsuite diversity beta-significance beta.tsv --metadata metadata.tsv --backend vegan --column site --method anosim2 --runtime docker -o anosim2.tsv
microsuite diversity_test --backend qiime2-adonis --distance-matrix distance.qza --metadata metadata.tsv --formula "batch + group / time_point" --permutations 999 --output adonis.qzv
microsuite functional_profile --backend picrust2 --table table.tsv --rep-seqs rep-seqs.fasta --output-dir functions --picrust2-database SC
microsuite report --backend native --run-dir runs/table-summary -o report.html
```

Metagenome assembly and binning examples:

```bash
microsuite assemble --backend megahit --read1 R1.fastq.gz --read2 R2.fastq.gz --output-dir assembly
microsuite assemble --backend mosh-megahit --reads reads.qza --output-contigs contigs.qza
microsuite assemble --backend idba-ud --reads reads.fa --output-dir assembly
microsuite bin --backend metabat2 --contigs contigs.fa --depth depth.tsv --output-dir bins
microsuite bin --backend mosh-metabat2 --contigs contigs.qza --alignment-maps reads-to-contigs-aln.qza --output-dir bins
microsuite bin --backend maxbin2 --contigs contigs.fa --abundance abundance.tsv --output-dir bins
microsuite bin --backend concoct --contigs contigs.fa --coverage coverage.tsv --output-dir bins
```

mothur backend examples (reference data is user-supplied; see [docs/mothur.md](mothur.md)):

```bash
microsuite cluster --backend mothur --rep-seqs contigs.fasta --reference-alignment silva.v4.align --count-table contigs.count_table --output-table table.tsv --output-rep-seqs rep-seqs.fasta --output-otu-list otu.list --output-count-table table.count_table --output-unique-seqs unique-seqs.fasta
microsuite tax_classify --backend mothur --rep-seqs unique-seqs.fasta --taxonomy-reference trainset.fasta --taxonomy-map trainset.tax --otu-list otu.list --count-table table.count_table -o taxonomy.tsv
microsuite workflow mothur --reads-dir reads/ --reference-alignment silva.v4.align --taxonomy-reference trainset.fasta --taxonomy-map trainset.tax --out runs/mothur-sop
```

Related QIIME 2 collections:

- MOSHPIT: https://moshpit.qiime2.org/en/stable/ for QIIME 2 whole-metagenome
  read-based and assembly-based workflows. The initial microsuite wrappers cover
  `mosh assembly assemble-megahit` and `mosh annotate bin-contigs-metabat`.

Conventions:

- command name = biological or analytical task
- `--backend` = tool or execution engine
- `--method` = statistical or transform method when needed
- `--run-dir` = optional runtime log directory for external-tool methods
- `--timeout` = optional external command timeout in seconds
- `--threads auto` = use detected CPU count minus one reserved core when supported
- external backends must fail clearly when executables are unavailable

The vegan backend delegates formula-based tests to `vegan::adonis2`. The
`anosim2` name is microsuite's stable compatibility name for vegan's
`vegan::anosim` function; vegan does not provide an upstream `anosim2` symbol.
For `adonis2`, `--formula` is the right-hand side only, and it can be combined
with `--strata` and `--permutations`. `--strata` accepts one metadata column or
a colon-separated interaction such as `batch:subject`; the latter is materialised
as a permutation block. This restricts permutations but does not fit an
lme4-style random effect, so `(1 | batch:subject)` is rejected explicitly.
The block is the joint interaction: rows shuffle only within the same batch
and subject. Unequal stratum sizes are supported, singleton strata cannot
move, and a design with only singleton strata is rejected. A warning identifies
coarse permutation spaces and formula columns with no within-stratum variation;
output records both requested and effective permutation counts. This interface
does not exchange whole subjects within batches; that requires an explicit
multi-level `permute::how()` design outside the current microsuite API.
Use `--runtime docker` with the `r-ecology` image when local R/vegan is not
installed.

The CLI is not responsible for large workflow scheduling. Use the Nextflow API
for full reproducible pipelines.

Runtime log directories contain `command.txt`, `stdout.log`, `stderr.log`,
`events.jsonl`, and `run.json`. These files are intended for debugging local
runs and for higher-level workflow systems that need lightweight provenance.

Deprecated compatibility commands may remain temporarily. For example,
`microsuite diffab ancombc` is kept as a compatibility alias for the preferred
`microsuite diff_abundance --backend ancombc` command.

For PICRUSt2 selectors, custom reference/trait contracts, outputs, coverage
semantics, Docker overrides, and migration from oldIMG, see
[PICRUSt2 functional profiling](picrust2.md).
