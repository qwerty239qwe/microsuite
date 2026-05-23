FROM ubuntu:24.04

LABEL org.opencontainers.image.title="mafft-fasttree"
LABEL org.opencontainers.image.description="Standalone MAFFT and FastTree phylogeny backend"

# Expected commands: mafft, FastTree
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates mafft fasttree && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["bash"]
