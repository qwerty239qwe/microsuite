# Toolbox model

`microsuite` is organized as a toolbox, not a library-first package.
It has three public APIs: Nextflow for full workflows, CLI for one-step tasks,
and Python SDK for programmatic table/statistics work.

The CLI is one public surface, alongside the Nextflow API and Python SDK:

- method-oriented commands such as `tax_classify --backend qiime2`
- `microsuite workflow ...` for end-to-end tasks
- `microsuite data ...` for demo datasets
- `microsuite qiime ...` for QIIME 2 artifact compatibility
- lower-level `import`, `diversity`, `ordination`, and `viz` commands for custom runs

The Python modules under `src/microsuite` are implementation details that
make the workflows testable and reusable. Stable notebook/script functions are
exposed through `microsuite.api`.

## Workflow convention

Each workflow writes a run directory containing:

- primary outputs such as `.h5ad`, `.tsv`, and figures
- `run.json` with inputs, parameters, outputs, and version
- `outputs.json` for quick programmatic discovery

Nextflow is the first-class API for full raw-read and external-tool workflows.
The Python CLI continues to own one-step local tasks, table import, summaries,
and visualization outputs.

## Method convention

Method-oriented commands are named by the biological or analytical task:

- `tax_classify --backend qiime2`
- `denoise --backend qiime2-dada2`
- `cluster --backend vsearch`
- `diversity_calc --backend qiime2`
- `normalize --backend native --method clr`
- `abundance --backend native --level genus`
- `shared_taxa --backend native --group body_site`
- `rarefy --backend native --depth 10000`
- `diff_abundance --backend ancombc`
- `report --backend native --run-dir runs/table-summary`
- future: `tax_classify --backend kraken2`
- future: `profile --backend metaphlan`

The `--backend` option selects the tool or engine. Backends that call external tools
must check for the external executable and fail with a clear setup message.
