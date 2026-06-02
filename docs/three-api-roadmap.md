# Three-API Roadmap

`microsuite` is a multi-environment microbiome toolbox. It exposes three
APIs for different users and execution contexts:

- Nextflow API: reproducible full workflows.
- CLI API: ergonomic one-step task commands.
- Python SDK: programmatic table/statistics functions.

Python remains installable because it powers the CLI and SDK, but the repository
is not limited to Python. R scripts, shell helpers, Dockerfiles, and Nextflow
workflows are first-class project assets.

## API Boundaries

### Nextflow API

Purpose: full reproducible pipelines.

Use for:

```text
FASTQ -> QC -> trim -> denoise/profile -> taxonomy -> phylogeny -> diversity -> report
```

Example:

```bash
nextflow run workflows/nextflow/main.nf \
  -profile docker \
  --workflow amplicon_qiime2 \
  --manifest manifest.tsv \
  --metadata metadata.tsv \
  --classifier classifier.qza \
  --outdir results
```

Owns:

- multi-step orchestration
- container/profile selection
- HPC/cloud execution
- resume and caching
- sample batching
- external tool environments

Does not own:

- native statistical implementations
- AnnData object logic
- interactive one-step commands

### CLI API

Purpose: one-step toolbox commands.

Examples:

```bash
microsuite qc --backend fastqc ...
microsuite denoise --backend qiime2-dada2 ...
microsuite tax_classify --backend qiime2 ...
microsuite normalize --backend native --method clr ...
microsuite diff_abundance --backend ancombc ...
```

Owns:

- task-oriented commands
- external tool wrappers
- local validation and error messages
- small manual runs
- workflow debugging

Does not own:

- large workflow scheduling
- container lifecycle beyond command checks
- HPC/cloud orchestration

### Python SDK

Purpose: programmatic downstream analysis.

Example:

```python
from microsuite.api import read_table, normalize_table, alpha_diversity

adata = read_table("table.h5ad")
adata = normalize_table(adata, method="clr")
alpha = alpha_diversity(adata, metric="shannon")
```

Owns:

- AnnData I/O
- native statistics
- normalization and transforms
- abundance tables
- diversity functions
- plotting helpers
- reusable logic for tests and notebooks

Does not own:

- QIIME2, DADA2, Kraken2, or R environment management
- Nextflow orchestration
- full external workflow execution

## Target Repository Layout

```text
.
src/microsuite/
  api/                  # public Python SDK facade
  cli/                  # CLI API
  methods/              # task backends used by CLI
  io/
  diversity/
  ordination/
  viz/
  workflows/            # Python-run mini workflows only
workflows/
  nextflow/
    main.nf
    nextflow.config
    modules/
    profiles/
containers/
  microsuite/
  qiime2-amplicon/
  r-diffab/
  kraken2/
scripts/
  r/
  shell/
examples/
docs/
tests/
```

## CLI Method Convention

Commands are named by biological or analytical task:

```bash
microsuite tax_classify --backend qiime2 ...
microsuite diversity_calc --backend qiime2 ...
microsuite normalize --backend native --method clr ...
```

Rules:

- `--backend` selects the tool or engine.
- `--method` is only used for a statistical or transform method.
- External backends check for the executable before running.
- Unsupported backends stay out of the public catalog until they have a working wrapper.
- Existing deprecated aliases can remain temporarily when needed.

## Task Surface

```text
qc
  - fastqc
  - multiqc
  - qiime2-demux

trim
  - fastp
  - cutadapt
  - trimmomatic
  - trim-galore
  - qiime2-cutadapt

denoise
  - qiime2-dada2
  - qiime2-deblur
  - dada2-r

cluster
  - vsearch
  - usearch
  - qiime2-vsearch

tax_classify
  - qiime2
  - kraken2
  - bracken
  - metaphlan
  - emu

phylogeny
  - qiime2
  - mafft-fasttree

normalize
  - native

abundance
  - native

shared_taxa
  - native

rarefy
  - native

diversity_calc
  - native
  - qiime2

beta_significance
  - qiime2
  - native

diff_abundance
  - ancombc
  - aldex2
  - maaslin2
  - lefse legacy

env_assoc
  - mantel
  - rda
  - cca
  - db-rda

network
  - native-correlation
  - sparcc
  - spieceasi
  - flashweave

functional_predict
  - picrust2
  - tax4fun2

functional_profile
  - humann

classify_samples
  - randomforest
  - xgboost optional

time_series
  - native

gamma_diversity
  - native

turnover
  - beta-turnover
  - taxa-turnover

visualize
  - native

report
  - native
```

## Implementation Phases

### Phase 1: Clean API Documentation and SDK Facade

Add:

```text
docs/api-cli.md
docs/api-python.md
docs/api-nextflow.md
docs/containers.md
src/microsuite/api/
```

Expose stable SDK functions:

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

Acceptance:

- SDK tests import from `microsuite.api`.
- Docs clearly distinguish Nextflow, CLI, and Python API.
- `README.md` and `docs/toolbox.md` are updated to match the three-API model.
- Existing CLI tests still pass.

### Phase 2: Differential Abundance Cleanup

Move the ANCOM-BC R code out of embedded Python strings.

Add:

```text
scripts/r/ancombc.R
src/microsuite/methods/diff_abundance.py
```

Preferred command:

```bash
microsuite diff_abundance --backend ancombc \
  --table table.h5ad \
  --group treatment \
  -o diff.tsv
```

Keep old compatibility command:

```bash
microsuite diffab ancombc ...
```

Acceptance:

- wrapper checks `Rscript`
- wrapper invokes packaged script path
- `ancombc`, `aldex2`, `maaslin2`, and `lefse` wrappers invoke packaged R scripts
- old command remains functional
- R implementation lives in `scripts/r/ancombc.R`, not in embedded Python strings

### Phase 3: Containers

Add skeletons:

```text
containers/microsuite/Dockerfile
containers/qiime2-amplicon/Dockerfile
containers/r-diffab/Dockerfile
containers/kraken2/Dockerfile
```

Acceptance:

- Dockerfiles include labels, purpose comments, and expected executables.
- Static tests verify files exist and mention expected tools.
- No container build required in the default unit test suite.

### Phase 4: Nextflow Skeleton

Add:

```text
workflows/nextflow/main.nf
workflows/nextflow/nextflow.config
workflows/nextflow/modules/
workflows/nextflow/profiles/
```

Initial workflow:

```text
amplicon_qiime2
```

Pipeline shape:

```text
manifest -> qiime2 dada2 -> taxonomy -> phylogeny -> diversity -> report
```

Acceptance:

- static tests verify core files and module names exist
- docs explain how to run with local/docker/singularity profiles
- full execution remains optional/manual until containers and sample data stabilize

### Phase 5: Reporting Layer

Add:

```bash
microsuite report --backend native \
  --run-dir runs/amplicon \
  -o report.html
```

Acceptance:

- report consumes `run.json` and optional `outputs.json`
- report renders a provenance summary of recorded inputs and outputs
- tests prove both metadata files are read from a tiny synthetic run directory

Deferred:

- summarize TSV contents and embed linked figures
- validate a formal Nextflow run-folder contract after the workflow emits
  `run.json`

## Review and TDD Process

For each major step:

1. Add failing tests first. For docs, containers, and Nextflow skeletons, use
   static existence/content tests before adding the files.
2. Implement the smallest code/docs needed.
3. Run:

```bash
uv run --extra dev pytest
uv run --extra dev ruff check .
uv run --extra dev ruff format --check .
uv run --extra dev ty check
```

4. Review docs to confirm the main goal is still represented.
5. Spawn a separate review agent to review the completed step.
6. Address review findings before moving to the next step.

Review agents should return findings with severity and file/line references.
For substantial phases, summarize accepted findings and fixes in
`docs/reviews/phase-N-review.md`.

If requirements are unclear or a design choice affects the public API, ask before
implementing.
