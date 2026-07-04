FROM rocker/r-ver:4.4.0

LABEL org.opencontainers.image.title="microsuite-dada2"
LABEL org.opencontainers.image.description="microsuite CLI + R DADA2 backend in one image (microsuite denoise --backend dada2-r)"

# Expected commands: microsuite, uv, Rscript
# Expected R packages: dada2

# System deps: Python to bootstrap uv, plus the compression/TLS libraries the
# DADA2 Bioconductor build links against.
RUN apt-get update && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*

# R DADA2 backend. Kept as its own layer so it caches independently of the
# (frequently changing) Python source below. dada2 is pinned by pinning the
# frozen Bioconductor release: 3.20 is the stable Bioconductor for the R 4.4
# base, so rebuilds resolve the same dada2 version. Bump with the R base image.
ARG BIOC_VERSION=3.20
RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    Rscript -e "BiocManager::install(version = '${BIOC_VERSION}', ask = FALSE, update = FALSE)" && \
    Rscript -e "BiocManager::install('dada2', ask = FALSE, update = FALSE)"

# microsuite CLI + SDK. `uv sync` installs the package so the bundled
# dada2_denoise.R ships inside the environment and denoise_dada2_r can find it.
RUN python3 -m pip install --no-cache-dir --break-system-packages uv
WORKDIR /opt/microsuite
COPY . /opt/microsuite
RUN uv sync

ENTRYPOINT ["uv", "run", "microsuite"]
