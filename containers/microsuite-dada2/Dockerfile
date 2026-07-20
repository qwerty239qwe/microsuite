FROM quay.io/biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2@sha256:64580af72a35a3ba07447c8b10e8bcf28a9f12ce7298e315c89df485d78354a2 AS dada2-runtime

FROM debian:bookworm-slim

# Biocontainers packages are installed into /usr/local. Copy that complete,
# mutually compatible R environment into a base that can also provide the
# process-inspection tools required by Nextflow.
COPY --from=dada2-runtime /usr/local/ /usr/local/

LABEL org.opencontainers.image.title="microsuite-dada2"
LABEL org.opencontainers.image.description="microsuite CLI + R DADA2 backend in one image (microsuite denoise --backend dada2-r)"

# Expected commands: microsuite, uv, Rscript
# Expected R packages: dada2

# DADA2 and its matching R/Bioconductor dependencies are supplied by the
# digest-pinned base image. Verify the expected package without resolving R
# repositories during this build.
RUN Rscript -e "library(dada2); stopifnot(as.character(packageVersion('dada2')) == '1.34.0')"

# microsuite CLI + SDK. Install uv via its standalone installer (base-agnostic;
# avoids the system-pip PEP-668 flag mismatch on the Ubuntu base). uv fetches its
# own Python for the venv, so `uv sync` installs the package and the packaged
# dada2_denoise.R that denoise_dada2_r invokes.
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"
WORKDIR /opt/microsuite
COPY . /opt/microsuite
RUN uv sync

ENTRYPOINT ["uv", "run", "microsuite"]
