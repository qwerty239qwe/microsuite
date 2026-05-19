# Containers

Container definitions are first-class toolbox assets because external microbiome
tools are difficult to install consistently.

Initial image skeletons:

```text
containers/microsuite/   Python CLI + SDK
containers/qiime2-amplicon/    QIIME 2 amplicon workflows
containers/r-diffab/           R differential-abundance tools
containers/kraken2/            Kraken2 taxonomy profiling
```

| Image | Purpose | Expected commands | Build status |
| --- | --- | --- | --- |
| `microsuite` | Python CLI and SDK runtime | `microsuite`, `uv` | skeleton |
| `qiime2-amplicon` | QIIME 2 amplicon backend | `qiime` | skeleton |
| `r-diffab` | R differential abundance backend | `Rscript`, `ANCOMBC` | skeleton |
| `kraken2` | Kraken2 taxonomy profiling | `kraken2`; planned: `bracken` | skeleton |

Default unit tests validate container files statically. GitHub Actions also
builds the lighter `microsuite` and `kraken2` images on every push and pull
request.

The CLI may check for external commands, but the Nextflow API should own
container/profile selection for full workflows.

Initial skeletons live under:

```text
containers/microsuite/Dockerfile
containers/qiime2-amplicon/Dockerfile
containers/r-diffab/Dockerfile
containers/kraken2/Dockerfile
```

Heavy images remain explicit validation steps because QIIME 2 and
R/Bioconductor images can be large and network-sensitive. Use the manual
GitHub Actions `build-heavy-containers=true` input to build `qiime2-amplicon`
and `r-diffab` in CI.

Build from the repository root so Dockerfiles that copy project files have the
right context:

```bash
docker build -f containers/microsuite/Dockerfile -t microsuite:local .
docker build -f containers/qiime2-amplicon/Dockerfile -t microsuite-qiime2-amplicon:local .
docker build -f containers/r-diffab/Dockerfile -t microsuite-r-diffab:local .
docker build -f containers/kraken2/Dockerfile -t microsuite-kraken2:local .
```
