FROM rocker/r-ver:4.4.0

LABEL org.opencontainers.image.title="r-dada2"
LABEL org.opencontainers.image.description="R DADA2 ASV inference backend"

# Expected commands: Rscript
# Expected R packages: dada2
#
# Pin dada2 by pinning the frozen Bioconductor release the image builds from:
# 3.20 is the stable Bioconductor for the R 4.4 base, so the installed dada2
# version is reproducible across rebuilds. Bump alongside the R base image.
ARG BIOC_VERSION=3.20
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install(version = '${BIOC_VERSION}', ask = FALSE, update = FALSE)" && \
    Rscript -e "BiocManager::install('dada2', ask = FALSE, update = FALSE)"

ENTRYPOINT ["Rscript"]
