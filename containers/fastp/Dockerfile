FROM ubuntu:24.04

LABEL org.opencontainers.image.title="fastp"
LABEL org.opencontainers.image.description="fastp trimming and read-filtering backend"

# Expected commands: fastp
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates fastp procps && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["fastp"]
