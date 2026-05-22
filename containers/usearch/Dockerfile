FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="usearch"
LABEL org.opencontainers.image.description="USEARCH 12 sequence search and clustering backend"

# Expected commands: usearch
RUN micromamba install -y -n base -c conda-forge -c bioconda usearch=12.0_beta && \
    micromamba clean -a -y

ENTRYPOINT ["usearch"]
