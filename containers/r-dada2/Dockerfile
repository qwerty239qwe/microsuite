FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="r-dada2"
LABEL org.opencontainers.image.description="R DADA2 ASV inference backend"

# Expected commands: Rscript
# Expected R packages: dada2
#
# The Bioconductor base provides an R/Bioconductor-compatible package stack and
# compile toolchain. Do not update its entire installed library while adding
# dada2; that makes builds much slower and less predictable. The trailing
# library() call turns a failed install into a hard build failure.
RUN Rscript -e "BiocManager::install('dada2', ask = FALSE, update = FALSE)" && \
    Rscript -e "library(dada2)"

ENTRYPOINT ["Rscript"]
