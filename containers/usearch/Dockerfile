FROM mambaorg/micromamba:2.0.5

LABEL org.opencontainers.image.title="usearch"
LABEL org.opencontainers.image.description="USEARCH 12 sequence search and clustering backend"

# Expected commands: usearch
USER root
RUN micromamba install -y -n base -c conda-forge -c bioconda usearch=12.0_beta && \
    printf '#!/bin/sh\nif [ -x /opt/conda/bin/usearch12 ]; then exec /opt/conda/bin/usearch12 "$@"; fi\nif [ -x /opt/conda/bin/usearch ]; then exec /opt/conda/bin/usearch "$@"; fi\necho "USEARCH executable not found in /opt/conda/bin" >&2\nexit 127\n' > /usr/local/bin/usearch && \
    chmod 755 /usr/local/bin/usearch && \
    micromamba clean -a -y

ENTRYPOINT ["/usr/local/bin/usearch"]
