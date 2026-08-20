# Tax4Fun2

microsuite 0.5.0 hardens Tax4Fun2 functional prediction around the upstream
Tax4Fun2 1.1.5 algorithm. It predicts functional potential from 16S nearest
reference genomes; it does not measure genes, transcripts, proteins, metabolites,
or activity.

## Reproducible execution

The recommended runtime is the pinned image. The large upstream reference data
is deliberately not embedded in the image and must be supplied separately.

```bash
microsuite functional_profile \
  --backend tax4fun2 \
  --table otu-table.tsv \
  --rep-seqs rep-seqs.fasta \
  --database Tax4Fun2_ReferenceData_v2 \
  --database-mode Ref99NR \
  --min-identity 0.97 \
  --output-dir tax4fun2-results \
  --runtime docker
```

The default image is
`ghcr.io/qwerty239qwe/microsuite/r-functional-tax4fun2:latest`. Override it with
`--image` or `MICROSUITE_R_FUNCTIONAL_TAX4FUN2_IMAGE`. Local execution requires
R, Tax4Fun2 exactly 1.1.5, jsonlite, `blastn`, and `makeblastdb` on `PATH`.

## Input contract

- The feature table is tab-delimited, with unique feature IDs in its first
  column and one or more uniquely named sample columns.
- Abundances must be numeric, finite, non-negative, and have a positive total
  in every sample.
- Representative sequences must be DNA FASTA records with unique identifiers.
- Table and FASTA feature-ID sets must match exactly. This prevents upstream
  merge behavior from silently discarding unpaired features.
- The selected `Ref99NR` or `Ref100NR` directory must contain its reference
  FASTA and compressed `.tbl.gz` profiles. `KEGG/ko.txt`, `KEGG/ko2ptw.txt`,
  and `KEGG/ptw.txt` are also required.

The reference directory is treated as input. BLAST indexes and all temporary
files are created under the staged output rather than beside the reference
FASTA.

## Output contract

| File | Meaning |
| --- | --- |
| `functional_prediction.tsv` | Relative KEGG ortholog prediction by sample. |
| `pathway_prediction.tsv` | Relative KEGG pathway prediction by sample. |
| `coverage.tsv` | Per-sample fractions of positive features and sequences represented by passing nearest-reference matches. |
| `tax4fun2_manifest.json` | Tax4Fun2 version, parameters, dimensions, matched-feature count, stable outputs, and MD5/size fingerprints for core and used reference files. |
| `upstream/` | BLAST logs/results, generated index, upstream logs, and original Tax4Fun2 text outputs. |
| `tax4fun2_container.json` | Docker image, engine, and digest when container execution is used. |

Treat low values in `coverage.tsv` as a warning that the predicted profile
represents only a limited part of the observed community. Always report the
Tax4Fun2 version, reference-data version, database mode, identity threshold,
copy-number normalization, pathway-normalization setting, and coverage.

## Scientific scope

Tax4Fun2 output is a reference-dependent inference from 16S sequences. Accuracy
depends on how well the sampled organisms are represented by sequenced and
functionally annotated genomes. It should be described as predicted functional
potential and should not be substituted for shotgun metagenomics,
metatranscriptomics, proteomics, or metabolomics when direct functional evidence
is required.
