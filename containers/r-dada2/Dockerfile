FROM rocker/r-ver:4.4.0

LABEL org.opencontainers.image.title="r-dada2"
LABEL org.opencontainers.image.description="R DADA2 ASV inference backend"

# Expected commands: Rscript
# Expected R packages: dada2
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install('dada2', ask = FALSE, update = FALSE)"

ENTRYPOINT ["Rscript"]
