FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-combatseq"
LABEL org.opencontainers.image.description="sva ComBat_seq batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: sva (ComBat_seq) (+ jsonlite for the params.json contract)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        bioconductor-sva && \
    micromamba clean -a -y

COPY src/microsuite/batch/r/combat_seq.R /opt/microsuite/combat_seq.R
COPY containers/r-batch-combatseq/smoke/ /opt/microsuite/smoke/

# Build-time smoke: a real correction, with an integer-output check. ComBat_seq
# returning floats would silently break every count-requiring downstream method.
RUN set -eux; \
    Rscript /opt/microsuite/combat_seq.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    Rscript -e "x <- read.delim('/tmp/microsuite-smoke-out.tsv', row.names=1); \
                stopifnot(all(abs(as.matrix(x) - round(as.matrix(x))) < 1e-8))"; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
