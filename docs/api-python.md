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
```

`read_table` and `write_table` currently support `.h5ad` files only. TSV, BIOM,
and QIIME 2 artifact import are available through the CLI and lower-level I/O
modules while the SDK surface stabilizes.

The SDK owns native table/statistics logic. It does not manage QIIME2, DADA2,
Kraken2, R, containers, or Nextflow execution environments.
