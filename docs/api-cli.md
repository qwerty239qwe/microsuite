# CLI API

The CLI API is for local, one-step microbiome tasks.

Use it when you want to run a single operation, debug a workflow step, or teach a
small example:

```bash
microsuite denoise --backend qiime2-dada2 --demux demux.qza --output-table table.qza --output-rep-seqs rep-seqs.qza --output-stats stats.qza --trunc-len 150
microsuite cluster --backend vsearch --table table.qza --rep-seqs rep-seqs.qza --output-table clustered-table.qza --output-rep-seqs clustered-rep-seqs.qza
microsuite normalize --backend native --method clr --table table.h5ad -o clr.h5ad
microsuite abundance --backend native --table table.h5ad --level genus -o abundance.tsv
microsuite tax_classify --backend qiime2 --rep-seqs rep-seqs.qza --classifier classifier.qza -o taxonomy.qza
microsuite diff_abundance --backend ancombc --table table.h5ad --group treatment -o diff.tsv
microsuite report --backend native --run-dir runs/table-summary -o report.html
```

Conventions:

- command name = biological or analytical task
- `--backend` = tool or execution engine
- `--method` = statistical or transform method when needed
- external backends must fail clearly when executables are unavailable

The CLI is not responsible for large workflow scheduling. Use the Nextflow API
for full reproducible pipelines.

Deprecated compatibility commands may remain temporarily. For example,
`microsuite diffab ancombc` is kept as a compatibility alias for the preferred
`microsuite diff_abundance --backend ancombc` command.
