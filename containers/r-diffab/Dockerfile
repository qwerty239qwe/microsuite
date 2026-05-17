FROM rocker/r-ver:4.4.0

LABEL org.opencontainers.image.title="r-diffab"
LABEL org.opencontainers.image.description="R differential abundance backend for ANCOMBC, ALDEx2, and related tools"

# Expected commands: Rscript
# Expected R packages: ANCOMBC
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install('ANCOMBC', ask = FALSE, update = FALSE)"

ENTRYPOINT ["Rscript"]
