FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-diffab"
LABEL org.opencontainers.image.description="R differential abundance backend for ANCOMBC, ALDEx2, MaAsLin2, and lefser"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R packages: ANCOMBC, ALDEx2, MaAsLin2, lefser
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        bioconductor-ancombc \
        bioconductor-aldex2 \
        bioconductor-maaslin2 \
        bioconductor-lefser && \
    micromamba clean -a -y && \
    Rscript -e "stopifnot(requireNamespace('ANCOMBC', quietly = TRUE)); stopifnot(requireNamespace('ALDEx2', quietly = TRUE)); stopifnot(requireNamespace('Maaslin2', quietly = TRUE)); stopifnot(requireNamespace('lefser', quietly = TRUE))"

COPY scripts/r/ancombc.R /opt/microsuite/ancombc.R
COPY scripts/r/aldex2.R /opt/microsuite/aldex2.R
COPY scripts/r/maaslin2.R /opt/microsuite/maaslin2.R
COPY scripts/r/lefse.R /opt/microsuite/lefse.R

ENTRYPOINT ["Rscript"]
