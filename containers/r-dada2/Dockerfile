FROM quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2@sha256:64580af72a35a3ba07447c8b10e8bcf28a9f12ce7298e315c89df485d78354a2

LABEL org.opencontainers.image.title="r-dada2"
LABEL org.opencontainers.image.description="R DADA2 ASV inference backend"

# Expected commands: Rscript
# Expected R packages: dada2
#
# DADA2 and its matching R/Bioconductor dependencies are already installed in
# this digest-pinned Biocontainers image. Verify the expected package eagerly;
# no R repository access is required while building this runtime wrapper.
RUN Rscript -e "library(dada2); stopifnot(as.character(packageVersion('dada2')) == '1.34.0')"

ENTRYPOINT ["Rscript"]
