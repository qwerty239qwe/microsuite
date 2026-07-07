FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="r-dada2"
LABEL org.opencontainers.image.description="R DADA2 ASV inference backend"

# Expected commands: Rscript
# Expected R packages: dada2
#
# The Bioconductor base pins R + Bioconductor 3.20 + a matching CRAN snapshot and
# ships the compile toolchain, so dada2 and its dependency tree install
# reproducibly. The trailing library() call turns a failed install (BiocManager
# only warns, so the layer would otherwise exit 0) into a hard build failure.
RUN Rscript -e "BiocManager::install('dada2', ask = FALSE, update = TRUE)" && \
    Rscript -e "library(dada2)"

ENTRYPOINT ["Rscript"]
