FROM bioconductor/bioconductor_docker:RELEASE_3_20

LABEL org.opencontainers.image.title="microsuite-dada2"
LABEL org.opencontainers.image.description="microsuite CLI + R DADA2 backend in one image (microsuite denoise --backend dada2-r)"

# Expected commands: microsuite, uv, Rscript
# Expected R packages: dada2

# R DADA2 backend. The Bioconductor base pins R + Bioconductor 3.20 + a matching
# CRAN snapshot and ships the compile toolchain, so dada2 installs reproducibly.
# The trailing library() call turns a failed install (BiocManager only warns)
# into a hard build failure instead of a silently broken image.
RUN Rscript -e "BiocManager::install('dada2', ask = FALSE, update = TRUE)" && \
    Rscript -e "library(dada2)"

# microsuite CLI + SDK. Install uv via its standalone installer (base-agnostic;
# avoids the system-pip PEP-668 flag mismatch on the Ubuntu base). uv fetches its
# own Python for the venv, so `uv sync` installs the package and the packaged
# dada2_denoise.R that denoise_dada2_r invokes.
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"
WORKDIR /opt/microsuite
COPY . /opt/microsuite
RUN uv sync

ENTRYPOINT ["uv", "run", "microsuite"]
