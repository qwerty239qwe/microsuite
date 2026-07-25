FROM condaforge/miniforge3:24.9.2-0

LABEL org.opencontainers.image.title="mothur"
LABEL org.opencontainers.image.description="mothur amplicon OTU clustering and classification backend"

# Expected commands: mothur
RUN conda install -y -c bioconda -c conda-forge mothur=1.48.5 && \
    conda clean -afy

ENTRYPOINT ["mothur"]
