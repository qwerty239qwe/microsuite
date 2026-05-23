FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="kraken2"
LABEL org.opencontainers.image.description="Kraken2 taxonomy profiling and Bracken abundance backend"

# Expected commands: kraken2, bracken
RUN micromamba install -y -n base -c conda-forge -c bioconda kraken2 bracken && \
    micromamba clean -a -y

ENTRYPOINT ["kraken2"]
