FROM bioconductor/bioconductor_docker:RELEASE_3_23

LABEL org.opencontainers.image.title="r-functional-tax4fun2"
LABEL org.opencontainers.image.description="Tax4Fun2 functional-prediction backend for microsuite"

ARG TAX4FUN2_COMMIT=91f08670021cf0921a7add3e67450795f23b92e5
ARG TAX4FUN2_SHA256=d74a69e3c5b2ca4aade42497bc65e0b5bc3cd55328a03dd2f20ff93a729a8a1f

# Tax4Fun2 1.1.5 is pinned by source commit. BLAST indexes are built in the
# writable output workspace at runtime, never beside the user reference data.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        ncbi-blast+ \
    && curl -L --fail --retry 3 \
        "https://github.com/fjossandon/Tax4Fun2/archive/${TAX4FUN2_COMMIT}.tar.gz" \
        -o /tmp/tax4fun2.tar.gz \
    && echo "${TAX4FUN2_SHA256}  /tmp/tax4fun2.tar.gz" | sha256sum -c - \
    && mkdir /tmp/tax4fun2 \
    && tar -xzf /tmp/tax4fun2.tar.gz --strip-components=1 -C /tmp/tax4fun2 \
    && R CMD INSTALL /tmp/tax4fun2 \
    && Rscript -e 'install.packages("jsonlite", repos="https://cloud.r-project.org"); stopifnot(as.character(packageVersion("Tax4Fun2")) == "1.1.5", as.character(packageVersion("jsonlite")) == "2.0.0")' \
    && rm -rf /tmp/tax4fun2 /tmp/tax4fun2.tar.gz /var/lib/apt/lists/*

COPY src/microsuite/functional/r/tax4fun2.R /opt/microsuite/tax4fun2.R
COPY containers/r-functional-tax4fun2/smoke/ /opt/microsuite/smoke/

# Build-time smoke: run real BLAST and Tax4Fun2 prediction against a minimal
# structurally valid reference, then verify planted KO direction and contracts.
RUN set -eux; \
    gzip /opt/microsuite/smoke/reference/Ref99NR/TESTA.tbl; \
    gzip /opt/microsuite/smoke/reference/Ref99NR/TESTB.tbl; \
    Rscript /opt/microsuite/tax4fun2.R \
        /opt/microsuite/smoke/rep-seqs.fasta \
        /opt/microsuite/smoke/table.tsv \
        /opt/microsuite/smoke/reference \
        1 Ref99NR 0.97 FALSE \
        /tmp/microsuite-tax4fun2-smoke; \
    test -s /tmp/microsuite-tax4fun2-smoke/functional_prediction.tsv; \
    test -s /tmp/microsuite-tax4fun2-smoke/pathway_prediction.tsv; \
    test -s /tmp/microsuite-tax4fun2-smoke/coverage.tsv; \
    test -s /tmp/microsuite-tax4fun2-smoke/tax4fun2_manifest.json; \
    awk -F '\t' 'NR == 1 { for (i=1;i<=NF;i++) { if ($i=="S_A") a=i; if ($i=="S_B") b=i } } $1=="K00001" { k1a=$a; k1b=$b } $1=="K00002" { k2a=$a; k2b=$b } END { exit !(a && b && k1a>k2a && k2b>k1b) }' \
        /tmp/microsuite-tax4fun2-smoke/functional_prediction.tsv; \
    awk -F '\t' 'NR > 1 { if ($2 != 1 || $3 != 1) exit 1; seen++ } END { exit !(seen == 2) }' \
        /tmp/microsuite-tax4fun2-smoke/coverage.tsv; \
    Rscript /opt/microsuite/tax4fun2.R \
        /opt/microsuite/smoke/rep-seqs.fasta \
        /opt/microsuite/smoke/table.tsv \
        /opt/microsuite/smoke/reference \
        /tmp/microsuite-tax4fun2-legacy \
        1 Ref99NR 0.97 FALSE; \
    test -s /tmp/microsuite-tax4fun2-legacy/functional_prediction.tsv; \
    ! Rscript /opt/microsuite/tax4fun2.R \
        /opt/microsuite/smoke/rep-seqs.fasta \
        /opt/microsuite/smoke/table.tsv \
        /opt/microsuite/smoke/reference \
        1 Ref99NR 0.97 MAYBE \
        /tmp/microsuite-tax4fun2-invalid; \
    rm -rf /tmp/microsuite-tax4fun2-smoke /tmp/microsuite-tax4fun2-legacy

ENTRYPOINT ["Rscript", "/opt/microsuite/tax4fun2.R"]
