# Moving Pictures QIIME 2 Parity

`microsuite workflow moving-pictures-qiime2` follows the major method sequence
from the QIIME 2 Moving Pictures tutorial, starting from raw EMP single-end
reads and training the Naive Bayes classifier.

Run:

```bash
microsuite workflow moving-pictures-qiime2 \
  --out runs/moving-pictures-qiime2 \
  --force
```

Requirements:

- activated QIIME 2 amplicon environment with `qiime` on `PATH`
- network access for `microsuite data fetch moving-pictures --full`
- enough local disk space for tutorial artifacts and visualizations

## Workflow Options

```bash
--out PATH
--force
--threads INT|auto
--timeout FLOAT
--qiime-command qiime
--include-deblur/--skip-deblur
--sampling-depth INT
--classifier-mode train
```

The default sampling depth is `1103`, matching the tutorial's common Moving
Pictures core-metrics setting. `--include-deblur` is off by default; Deblur is
available as a method wrapper but the parity workflow follows the DADA2 path.

## Tutorial Mapping

| Tutorial method | microsuite command or workflow step | Expected output |
| --- | --- | --- |
| Metadata tabulation | `microsuite metadata_tabulate --backend qiime2` | `sample-metadata.qzv` |
| EMP import | `microsuite qiime_import --backend qiime2-emp-single-end` | `emp-single-end-sequences.qza` |
| EMP demultiplexing | `microsuite demux --backend qiime2-emp-single` | `demux.qza`, `demux-details.qza` |
| Demux summary | `microsuite qc --backend qiime2-demux` | `demux.qzv` |
| DADA2 single-end | `microsuite denoise --backend qiime2-dada2 --trim-left 0 --trunc-len 120` | `table.qza`, `rep-seqs.qza`, `denoising-stats.qza` |
| Denoising stats tabulation | `microsuite metadata_tabulate --backend qiime2` | `denoising-stats.qzv` |
| Feature table summary | `microsuite feature_summarize --backend qiime2 --mode summarize` | `table.qzv` |
| Representative sequences | `microsuite feature_summarize --backend qiime2 --mode tabulate-seqs` | `rep-seqs.qzv` |
| MAFFT/FastTree phylogeny | `microsuite phylogeny --backend qiime2-mafft-fasttree` | aligned sequences, masked alignment, unrooted tree, rooted tree |
| Core diversity | `microsuite diversity_core --backend qiime2-core-metrics-phylogenetic` | QIIME core-metrics directory |
| Alpha significance | `microsuite diversity_test --backend qiime2-alpha-group-significance` | alpha group-significance `.qzv` |
| Beta significance | `microsuite diversity_test --backend qiime2-beta-group-significance` | beta group-significance `.qzv` |
| Emperor plot | `microsuite ordination_plot --backend qiime2-emperor` | Emperor `.qzv` |
| Alpha rarefaction | `microsuite rarefaction --backend qiime2-alpha-rarefaction` | alpha-rarefaction `.qzv` |
| Classifier training | `microsuite tax_train --backend qiime2-naive-bayes` | trained classifier `.qza` |
| Taxonomy classification | `microsuite tax_classify --backend qiime2` | taxonomy `.qza` |
| Taxonomy tabulation | `microsuite metadata_tabulate --backend qiime2` | taxonomy `.qzv` |
| Taxa barplot | `microsuite tax_barplot --backend qiime2` | taxa barplot `.qzv` |
| Gut sample filtering | `microsuite feature_filter --backend qiime2-filter-samples` | gut-only table `.qza` |
| ANCOM-BC | `microsuite diff_abundance --backend qiime2-ancombc` | differentials `.qza` |
| DA barplot | `microsuite diff_viz --backend qiime2-da-barplot` | DA barplot `.qzv` |
| Level-6 collapse | `microsuite tax_collapse --backend qiime2 --level 6` | level-6 table `.qza` |
| Level-6 ANCOM-BC and plot | `diff_abundance` then `diff_viz` | level-6 differentials and `.qzv` |

## Outputs

The workflow writes:

- `data/`: downloaded tutorial inputs
- `qiime2/`: QIIME 2 artifacts
- `visualizations/`: QIIME 2 visualizations
- `core-metrics-results/`: diversity outputs
- `runtime/NN-step-name/`: `command.txt`, `stdout.log`, `stderr.log`,
  `events.jsonl`, and `run.json` for each QIIME command
- `run.json`: ordered workflow provenance
- `report.html`: native microsuite summary of executed steps

## Known Limitations

- Default tests mock QIIME 2 command execution; real QIIME 2 integration is
  manual or environment-gated.
- The workflow is QIIME-wrapper-first for 0.1.0 and does not reimplement DADA2,
  MAFFT/FastTree, diversity, taxonomy, or ANCOM-BC natively.
- Tutorial URLs are pinned to the QIIME 2 2024.10 tutorial data family until a
  release-pinned container/runtime is finalized.
- Quick precomputed-artifact demos remain separate:

```bash
microsuite workflow moving-pictures --out runs/moving-pictures --force
microsuite data fetch moving-pictures -o data/moving-pictures-real --full
```
