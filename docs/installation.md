# Installation

`microsuite` has a small Python core and many optional external backends. Pick
the setup that matches your work. You do not need every tool unless you plan to
run every backend.

## Which Setup Do I Need?

| User type | Install | Best for |
| --- | --- | --- |
| Python SDK users | Python 3.11/3.12 + optional Python extras | Notebooks, scripts, AnnData table analysis |
| CLI users | Python core plus only the external tools you call | One-step local commands |
| Docker users | Docker images | Avoiding local installs of heavy bioinformatics tools |
| Nextflow and HPC users | Nextflow plus local/Docker/Singularity profile | Full reproducible pipelines |
| Cloud/Spot VM users | Docker image plus cloud orchestration scripts | Large BioProject-scale jobs |
| Developers | `uv sync --extra dev` | Tests, linting, editing the repo |

## Developer/Local Install

Use this when working from a clone of the repository:

```bash
uv sync --extra dev
uv run microsuite --help
uv run pytest
```

This installs the Python package, test tools, and development tooling. It does
not install QIIME 2, R packages, Kraken2 databases, SRA Toolkit, or other
external bioinformatics runtimes.

## Python SDK Users

Use the SDK when you already have feature tables and want notebook/script access
to native Python analysis.

```bash
uv sync
uv run python -c "from microsuite.api import alpha_diversity; print(alpha_diversity)"
```

Optional Python extras:

```bash
uv sync --extra biom
uv sync --extra qza
uv sync --extra all
```

What those extras mean:

| Extra | Adds | Use when |
| --- | --- | --- |
| `biom` | `biom-format` | Reading BIOM feature tables |
| `qza` | `biom-format` | Reading QIIME 2 artifacts containing BIOM tables |
| `all` | All Python extras | You want all Python-side file compatibility |
| `dev` | pytest, ruff, ty | You are developing or testing the repo |

Minimal SDK smoke test:

```bash
uv run python - <<'PY'
from microsuite.api import read_table
print("microsuite SDK import works")
PY
```

## CLI Users

Use this path for local one-step commands.

```bash
uv sync
uv run microsuite --help
uv run microsuite workflow list
```

Native table/statistics commands work from the Python environment. Commands
that wrap external tools need those tools installed on `PATH`.

Common examples:

| Command family | Also install |
| --- | --- |
| `qc --backend fastqc` | FastQC |
| `qc --backend multiqc` | MultiQC |
| `trim --backend fastp` | fastp |
| `trim --backend cutadapt` | Cutadapt |
| `cluster --backend vsearch` | VSEARCH |
| `denoise --backend dada2-r` | R, `Rscript`, DADA2 |
| `diversity alpha --metric breakaway` | R, `Rscript`, Breakaway |
| `diversity alpha --metric inext` | R, `Rscript`, iNEXT |
| QIIME 2 backends | Activated QIIME 2 environment with required plugins |
| Kraken2/Bracken | Tool binaries plus prepared databases |

If an external dependency is missing, the command should fail with an actionable
message naming the missing runtime or package.

## Docker Users

Use Docker when you do not want to install every bioinformatics tool locally.

Build from the repository root:

```bash
docker build -f containers/microsuite/Dockerfile -t microsuite:local .
docker run --rm microsuite:local --help
```

Tool-specific containers live under `containers/` and are documented in
[containers.md](containers.md). Build only the images you need:

```bash
docker build -f containers/fastqc/Dockerfile -t microsuite-fastqc:local .
docker build -f containers/vsearch/Dockerfile -t microsuite-vsearch:local .
docker build -f containers/r-diffab/Dockerfile -t microsuite-r-diffab:local .
```

Mount your data into the container:

```bash
docker run --rm \
  -v "$PWD/data:/data" \
  -v "$PWD/results:/results" \
  microsuite:local \
  microsuite --help
```

## Nextflow And HPC Users

Use Nextflow for complete, resumable workflows.

Requirements depend on profile:

| Profile | Requirements |
| --- | --- |
| `local` | Tools such as `fastqc`, `multiqc`, and `qiime` available on `PATH` |
| `docker` | Docker available on the machine or executor |
| `singularity` | Singularity/Apptainer plus matching image files |

Smoke test the workflow entry point:

```bash
nextflow run workflows/nextflow/main.nf -profile docker --help
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

See [api-nextflow.md](api-nextflow.md) for the manifest contract and workflow
parameters.

> Running many samples? See [multisample runs and concurrency](multisample.md).

## R And QIIME 2 Backends

Python extras do not install R packages or QIIME 2.

For R-backed methods, verify:

```bash
Rscript -e 'stopifnot(requireNamespace("breakaway", quietly = TRUE))'
Rscript -e 'stopifnot(requireNamespace("iNEXT", quietly = TRUE))'
```

For QIIME 2-backed methods, activate your QIIME 2 environment first:

```bash
qiime --help
uv run microsuite tax_classify --help
```

Keep database paths explicit. Tools such as Kraken2, Bracken, MetaPhlAn, EMU,
PICRUSt2, HUMAnN, and Tax4Fun2 generally need large user-supplied databases.

## Cloud/Spot VM Users

Use this for large jobs where cheap interrupted compute is acceptable.

These scripts are optional manual validation infrastructure. They create
billable cloud resources and are not run by default CI.

The PRJNA321534 runner has its own cloud kit:

```text
cloud/prjna321534/
cloud/gcp/
containers/prjna321534-alpha/
```

Typical flow:

```bash
# 1. Build/push the image with GitHub Actions:
#    PRJNA321534 Cloud Image

# 2. Start a small pilot Spot VM:
export PROJECT_ID="my-gcp-project"
export BUCKET="microsuite-prjna321534-$USER"
export IMAGE="ghcr.io/qwerty239qwe/microsuite-prjna321534-alpha:latest"
export MAX_RUNS="10"
bash cloud/gcp/create_spot_vm.sh

# 3. Download completed results:
bash cloud/gcp/download_results.sh

# 4. Delete the VM:
bash cloud/gcp/delete_spot_vm.sh
```

Spot VMs can be interrupted. Keep work and results in object storage and run a
small pilot before launching a full BioProject run.

## Verification Checklist

Run the checks that match your setup:

```bash
uv run microsuite --help
uv run microsuite workflow list
uv run pytest
uv run ruff check .
uv run ty check
docker run --rm microsuite:local --help
nextflow run workflows/nextflow/main.nf -profile docker --help
```

For external backends, also verify the external command itself:

```bash
fastqc --version
vsearch --version
Rscript --version
qiime --help
```

## Troubleshooting

| Symptom | Likely cause | Fix |
| --- | --- | --- |
| `ModuleNotFoundError: biom` | Missing Python extra | Run `uv sync --extra biom` or `uv sync --extra qza` |
| `Rscript` not found | R is not installed or not on `PATH` | Install R or use a container |
| Breakaway/iNEXT errors | R package missing | Install `breakaway` or `iNEXT` in the R library used by `Rscript` |
| QIIME 2 command fails | QIIME 2 environment not active | Activate QIIME 2 before running the wrapper |
| Kraken2/Bracken runs but no useful output | Missing or mismatched database | Provide a prepared database matching the command |
| Nextflow Docker run cannot start | Docker daemon/profile issue | Test `docker run hello-world` and check the selected profile |
