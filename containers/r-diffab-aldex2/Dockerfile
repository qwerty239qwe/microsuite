FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-diffab-aldex2"
LABEL org.opencontainers.image.description="ALDEx2 differential-abundance R backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: ALDEx2
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        bioconductor-aldex2 && \
    micromamba clean -a -y

COPY src/microsuite/diffab/r/aldex2.R /opt/microsuite/aldex2.R
COPY containers/r-diffab-aldex2/smoke/ /opt/microsuite/smoke/

# Build-time smoke: run a tiny real dataset through the backend (not just a
# package import), print diagnostics, and fail the build if it can't produce a
# non-empty result.
RUN set -eux; \
    Rscript /opt/microsuite/aldex2.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        group \
        /opt/microsuite/smoke_out.tsv; \
    test -s /opt/microsuite/smoke_out.tsv; \
    rm -f /opt/microsuite/smoke_out.tsv

ENTRYPOINT ["Rscript"]
