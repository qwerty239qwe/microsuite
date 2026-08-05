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
# build unless it emits a table with every input feature and sample.
RUN set -eux; \
    Rscript /opt/microsuite/mmuphin.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
