# Data model

The canonical object is `AnnData`.

- `adata.X`: sample x feature count matrix
- `adata.obs`: sample metadata indexed by sample ID
- `adata.var`: feature metadata indexed by feature ID
- `adata.var["taxonomy"]`: raw taxonomy string when provided
- `adata.uns["microsuite"]`: importer provenance

Taxonomy strings are parsed into these columns when possible:

- `kingdom`
- `phylum`
- `class`
- `order`
- `family`
- `genus`
- `species`
