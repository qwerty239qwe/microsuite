FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-diffab-ancombc"
LABEL org.opencontainers.image.description="ANCOM-BC differential-abundance R backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: ANCOMBC (+ jsonlite for the params.json contract)
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        bioconductor-ancombc && \
    micromamba clean -a -y

COPY src/microsuite/diffab/r/ancombc.R /opt/microsuite/ancombc.R
COPY containers/r-diffab-ancombc/smoke/ /opt/microsuite/smoke/

# Build-time smoke: run a tiny real dataset through the backend (not just a
# package import), print diagnostics, and fail the build if it can't produce a
# non-empty result.
RUN set -eux; \
    Rscript /opt/microsuite/ancombc.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /opt/microsuite/smoke_out.tsv; \
    test -s /opt/microsuite/smoke_out.tsv; \
    rm -f /opt/microsuite/smoke_out.tsv

ENTRYPOINT ["Rscript"]
