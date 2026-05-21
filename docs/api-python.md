# Python SDK API

The Python SDK is for notebooks, custom scripts, and reusable downstream
analysis. It works on AnnData tables and exposes stable function names from
`microsuite.api`.

Example:

```python
from microsuite.api import read_table, normalize_table, abundance_table

adata = read_table("table.h5ad")
relative = normalize_table(adata, method="relative")
genus = abundance_table(relative, level="genus")
```

Initial SDK surface:

```python
read_table(path)
write_table(adata, path)
normalize_table(adata, method="relative")
abundance_table(adata, level="genus")
shared_taxa_table(adata, level="genus", group="body_site")
rarefy_table(adata, depth=10000)
alpha_diversity(adata, metric="shannon")
beta_diversity(adata, metric="bray-curtis")
pcoa(distance_matrix)
qc(backend="fastqc", inputs=[path], output_dir=qc_dir)
trim(backend="cutadapt", read1=read, output1=trimmed, adapter=adapter)
trim(backend="trimmomatic", read1=read, output1=trimmed, trimmomatic_steps=steps)
trim(backend="trim-galore", read1=read, output1=trimmed, trim_galore_version="auto")
denoise(backend="dada2-r", demux=reads_dir, output_table=table, output_rep_seqs=rep_seqs, output_stats=stats)
qc_filter(backend="qiime2-filter-reads", demux=demux, database=index, output=filtered)
decontam(backend="qiime2-decontam", table=table, metadata=metadata, output=scores)
evaluate(backend="qiime2-taxonomy", expected_taxa=expected, observed_taxa=observed, output=viz)
```

`read_table` and `write_table` currently support `.h5ad` files only. TSV, BIOM,
and QIIME 2 artifact import are available through the CLI and lower-level I/O
modules while the SDK surface stabilizes.

The SDK owns native table/statistics logic and can launch selected external
tool methods, such as FastQC, through the same backend names used by the CLI.
It does not install or activate QIIME 2, DADA2, Kraken2, R, containers, or
Nextflow execution environments; those tools must already be available in the
runtime environment.

External-tool methods accept optional `run_dir` and `timeout` arguments. When
`run_dir` is supplied, microsuite writes `command.txt`, `stdout.log`,
`stderr.log`, `events.jsonl`, and `run.json` for debugging and provenance.
Threaded methods accept a positive integer, and the method-oriented wrappers
also accept `threads="auto"` where the backend supports threading.

FastQC example:

```python
from pathlib import Path

from microsuite.api import qc

qc(
    backend="fastqc",
    inputs=[Path("sample.fastq.gz")],
    output_dir=Path("qc"),
    threads=4,
    extract=True,
    run_dir=Path("runs/fastqc/sample"),
)
```

Cutadapt example:

```python
from pathlib import Path

from microsuite.api import trim

trim(
    backend="cutadapt",
    read1=Path("sample_R1.fastq.gz"),
    output1=Path("trimmed_R1.fastq.gz"),
    adapter="AGATCGGAAGAGC",
    quality_cutoff="20",
    minimum_length="100",
    json_report=Path("cutadapt.json"),
    threads=4,
)
```

Trimmomatic example:

```python
from pathlib import Path

from microsuite.api import trim

trim(
    backend="trimmomatic",
    read1=Path("sample_R1.fastq.gz"),
    output1=Path("trimmed_R1.fastq.gz"),
    trimmomatic_steps=["SLIDINGWINDOW:4:20", "MINLEN:100"],
    threads=4,
)
```

R/DADA2 example:

```python
from pathlib import Path

from microsuite.api import denoise

denoise(
    backend="dada2-r",
    demux=Path("trimmed-fastq"),
    output_table=Path("table.tsv"),
    output_rep_seqs=Path("rep-seqs.fasta"),
    output_stats=Path("stats.tsv"),
    paired=True,
    trunc_len_f=151,
    trunc_len_r=149,
    threads=4,
)
```

QIIME 2 quality-control example:

```python
from pathlib import Path

from microsuite.api import qc_filter

qc_filter(
    backend="qiime2-filter-reads",
    demux=Path("demux.qza"),
    database=Path("human-bowtie2-index.qza"),
    output=Path("filtered.qza"),
    threads=8,
)
```
