FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-mmuphin"
LABEL org.opencontainers.image.description="MMUPHin adjust_batch batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: MMUPHin (+ jsonlite for the params.json contract)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        bioconductor-mmuphin && \
    micromamba clean -a -y

COPY src/microsuite/batch/r/mmuphin.R /opt/microsuite/mmuphin.R
COPY containers/r-batch-mmuphin/smoke/ /opt/microsuite/smoke/

# Build-time smoke: run a real correction, not a package import, and fail the
# build unless it emits a table with every input feature and sample. Also
# assert the declared `relative` scale (backends.py): mmuphin.R TSS-
# normalizes adjust_batch()'s count-scale output before writing, so every
# sample column must sum to ~1.
RUN set -eux; \
    Rscript /opt/microsuite/mmuphin.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    Rscript -e "x <- read.delim('/tmp/microsuite-smoke-out.tsv', row.names=1); \
                stopifnot(all(abs(colSums(as.matrix(x)) - 1) < 1e-6))"; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
