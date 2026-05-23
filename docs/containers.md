# Containers

Container definitions are first-class toolbox assets because external microbiome
tools are difficult to install consistently.

Initial image skeletons:

```text
containers/microsuite/   Python CLI + SDK
containers/fastqc/            FastQC raw-read quality reports
containers/multiqc/           MultiQC aggregate quality reports
containers/fastp/             fastp trimming and read filtering
containers/cutadapt/          Cutadapt adapter and primer trimming
containers/trimmomatic/       Trimmomatic read trimming
containers/trim-galore/       Trim Galore wrapper around Cutadapt
containers/vsearch/           VSEARCH sequence clustering
containers/usearch/           USEARCH 12 sequence search and clustering
containers/mafft-fasttree/     MAFFT/FastTree phylogeny
containers/metaphlan/          MetaPhlAn marker-gene profiling
containers/qiime2-amplicon/    QIIME 2 amplicon workflows
containers/r-dada2/            R DADA2 ASV inference
containers/r-diffab/           R differential-abundance tools
containers/kraken2/            Kraken2 taxonomy profiling and Bracken abundance
```

| Image | Purpose | Expected commands | Build status |
| --- | --- | --- | --- |
| `microsuite` | Python CLI and SDK runtime | `microsuite`, `uv` | skeleton |
| `fastqc` | Raw-read quality reports | `fastqc` | implemented |
| `multiqc` | Aggregate quality reports | `multiqc` | implemented |
| `fastp` | Adapter trimming, quality filtering, and reports | `fastp` | implemented |
| `cutadapt` | Adapter and primer trimming | `cutadapt` | implemented |
| `trimmomatic` | Sliding-window, length, quality, and adapter trimming | `trimmomatic` | implemented |
| `trim-galore` | Trim Galore trimming wrapper | `trim_galore`, `cutadapt` | implemented |
| `vsearch` | VSEARCH sequence clustering | `vsearch` | implemented |
| `usearch` | USEARCH 12 sequence search and clustering | `usearch` | implemented |
| `mafft-fasttree` | Standalone alignment and phylogeny | `mafft`, `FastTree` | implemented |
| `metaphlan` | MetaPhlAn marker-gene profiling | `metaphlan` | implemented |
| `qiime2-amplicon` | QIIME 2 amplicon backend | `qiime` | skeleton |
| `r-dada2` | R DADA2 ASV inference | `Rscript`, `dada2` | implemented |
| `r-diffab` | R differential abundance backend | `Rscript`, `ANCOMBC` | skeleton |
| `kraken2` | Kraken2 taxonomy profiling and Bracken abundance | `kraken2`, `bracken` | implemented |

Default unit tests validate container files statically. The separate Docker
GitHub Actions workflow builds the lighter `microsuite`, `fastqc`, trimming,
MultiQC, VSEARCH, USEARCH, MAFFT/FastTree, and `kraken2` images when container
files change or when it is run manually. Lightweight images run tiny
mounted-input examples after they build.

The CLI may check for external commands, but the Nextflow API should own
container/profile selection for full workflows.

The `fastqc` image is the first concrete external-tool runtime. It installs
FastQC 0.12.1 with Java and the upstream `fastqc` wrapper, and CI smoke-tests
`fastqc --version`. The current Nextflow `FASTQC` module is still a workflow
placeholder until the raw-read manifest contract is finalized.

Initial skeletons live under:

```text
containers/microsuite/Dockerfile
containers/fastqc/Dockerfile
containers/multiqc/Dockerfile
containers/fastp/Dockerfile
containers/cutadapt/Dockerfile
containers/trimmomatic/Dockerfile
containers/trim-galore/Dockerfile
containers/vsearch/Dockerfile
containers/usearch/Dockerfile
containers/mafft-fasttree/Dockerfile
containers/metaphlan/Dockerfile
containers/qiime2-amplicon/Dockerfile
containers/r-dada2/Dockerfile
containers/r-diffab/Dockerfile
containers/kraken2/Dockerfile
```

Heavy images remain explicit validation steps because QIIME 2 and
R/Bioconductor images can be large and network-sensitive. Use the manual
GitHub Actions `build-heavy-containers=true` input to build `metaphlan`,
`qiime2-amplicon`, `r-dada2`, and `r-diffab` in CI.

Build from the repository root so Dockerfiles that copy project files have the
right context:

```bash
docker build -f containers/microsuite/Dockerfile -t microsuite:local .
docker build -f containers/fastqc/Dockerfile -t microsuite-fastqc:local .
docker build -f containers/multiqc/Dockerfile -t microsuite-multiqc:local .
docker build -f containers/fastp/Dockerfile -t microsuite-fastp:local .
docker build -f containers/cutadapt/Dockerfile -t microsuite-cutadapt:local .
docker build -f containers/trimmomatic/Dockerfile -t microsuite-trimmomatic:local .
docker build -f containers/trim-galore/Dockerfile -t microsuite-trim-galore:local .
docker build -f containers/vsearch/Dockerfile -t microsuite-vsearch:local .
docker build -f containers/usearch/Dockerfile -t microsuite-usearch:local .
docker build -f containers/mafft-fasttree/Dockerfile -t microsuite-mafft-fasttree:local .
docker build -f containers/metaphlan/Dockerfile -t microsuite-metaphlan:local .
docker build -f containers/qiime2-amplicon/Dockerfile -t microsuite-qiime2-amplicon:local .
docker build -f containers/r-dada2/Dockerfile -t microsuite-r-dada2:local .
docker build -f containers/r-diffab/Dockerfile -t microsuite-r-diffab:local .
docker build -f containers/kraken2/Dockerfile -t microsuite-kraken2:local .
```
