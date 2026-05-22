FROM ubuntu:24.04

LABEL org.opencontainers.image.title="trimmomatic"
LABEL org.opencontainers.image.description="Trimmomatic read trimming backend"

# Expected commands: trimmomatic, java, openjdk
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates trimmomatic && \
    printf '#!/bin/sh\nexec java -jar /usr/share/java/trimmomatic.jar "$@"\n' \
      > /usr/local/bin/trimmomatic && \
    chmod 755 /usr/local/bin/trimmomatic && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["trimmomatic"]
