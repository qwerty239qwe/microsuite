FROM ubuntu:24.04

LABEL org.opencontainers.image.title="vsearch"
LABEL org.opencontainers.image.description="VSEARCH sequence clustering backend"

# Expected commands: vsearch
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates vsearch && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["vsearch"]
