FROM ubuntu:24.04

LABEL org.opencontainers.image.title="trim-galore"
LABEL org.opencontainers.image.description="Trim Galore adapter and quality trimming backend"

# Expected commands: trim_galore, cutadapt
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates trim-galore && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["trim_galore"]
