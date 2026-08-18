FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL org.opencontainers.image.title="r-diffab-maaslin3"
LABEL org.opencontainers.image.description="MaAsLin 3 differential-abundance R backend for microsuite"

# Expected commands: Rscript
# Expected R package: maaslin3 1.4.0
RUN Rscript -e 'BiocManager::install(c("maaslin3", "jsonlite"), ask = FALSE, update = FALSE); stopifnot(as.character(packageVersion("maaslin3")) == "1.4.0")'

COPY src/microsuite/diffab/r/maaslin3.R /opt/microsuite/maaslin3.R
COPY containers/r-diffab-maaslin3/smoke/ /opt/microsuite/smoke/

# Build-time smoke: fit both model families and require distinct result tables.
RUN set -eux; \
    Rscript /opt/microsuite/maaslin3.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-maaslin3-smoke; \
    test -s /tmp/microsuite-maaslin3-smoke/abundance_results.tsv; \
    test -s /tmp/microsuite-maaslin3-smoke/prevalence_results.tsv; \
    awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "model") model = i; next } { if ($model != "abundance") bad = 1; seen = 1 } END { exit !(model && seen && !bad) }' /tmp/microsuite-maaslin3-smoke/abundance_results.tsv; \
    awk -F '\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "model") model = i; next } { if ($model != "prevalence") bad = 1; seen = 1 } END { exit !(model && seen && !bad) }' /tmp/microsuite-maaslin3-smoke/prevalence_results.tsv; \
    rm -rf /tmp/microsuite-maaslin3-smoke

ENTRYPOINT ["Rscript"]
