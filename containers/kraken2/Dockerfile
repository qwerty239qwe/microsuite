FROM ubuntu:24.04

LABEL org.opencontainers.image.title="kraken2"
LABEL org.opencontainers.image.description="Kraken2 taxonomy profiling and Bracken abundance backend"

# Expected commands: kraken2, bracken
RUN apt-get update && \
    apt-get install -y --no-install-recommends kraken2 bracken ca-certificates && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["kraken2"]
