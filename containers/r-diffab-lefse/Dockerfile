FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL org.opencontainers.image.title="r-diffab-lefse"
LABEL org.opencontainers.image.description="LEfSe (lefser) differential-abundance R backend for microsuite"

# Expected commands: Rscript
# Expected R package: lefser 1.22.0
RUN Rscript -e 'BiocManager::install(c("lefser", "jsonlite"), ask = FALSE, update = FALSE); stopifnot(as.character(packageVersion("lefser")) == "1.22.0")'

COPY src/microsuite/diffab/r/lefse.R /opt/microsuite/lefse.R
COPY containers/r-diffab-lefse/smoke/ /opt/microsuite/smoke/

# Build-time smoke: run a tiny real dataset through the backend (not just a
# package import), print diagnostics, and fail the build if it can't produce a
# non-empty result.
RUN set -eux; \
    Rscript /opt/microsuite/lefse.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    awk -F '\t' 'NR == 1 { exit !(NF == 2 && $1 == "features" && $2 == "scores") }' /tmp/microsuite-smoke-out.tsv; \
    grep -q '^F0\|^F1\|^F2' /tmp/microsuite-smoke-out.tsv; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
