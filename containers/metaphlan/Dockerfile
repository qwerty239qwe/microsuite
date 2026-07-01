FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="metaphlan"
LABEL org.opencontainers.image.description="MetaPhlAn marker-gene taxonomic profiling backend"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: metaphlan
RUN micromamba install -y -n base -c conda-forge -c bioconda metaphlan && \
    micromamba clean -a -y

ENTRYPOINT ["metaphlan"]
