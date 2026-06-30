FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="microsuite-picrust2"
LABEL org.opencontainers.image.description="microsuite CLI + PICRUSt2 for functional profiling of 16S OTUs"

ENV PATH="/opt/conda/bin:${PATH}"

# PICRUSt2 bundles its own reference data; install it alongside Python so the
# microsuite CLI and picrust2_pipeline.py share one environment.
# Expected commands: picrust2_pipeline.py, microsuite
RUN micromamba install -y -n base -c conda-forge -c bioconda picrust2 pip && \
    micromamba clean -a -y

USER root
WORKDIR /opt/microsuite
COPY . /opt/microsuite
RUN pip install --no-cache-dir ".[biom,qza]"

ENTRYPOINT ["microsuite"]
