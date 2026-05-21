FROM ubuntu:24.04

ARG FASTQC_VERSION=0.12.1
ARG FASTQC_SHA256=5f4dba8780231a25a6b8e11ab2c238601920c9704caa5458d9de559575d58aa7

LABEL org.opencontainers.image.title="fastqc"
LABEL org.opencontainers.image.description="FastQC quality-report backend for raw sequencing reads"

# Expected commands: fastqc, java, perl, openjdk
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      curl \
      libfontconfig1 \
      libfreetype6 \
      libharfbuzz0b \
      openjdk-21-jre-headless \
      perl \
      unzip && \
    curl -fsSL \
      "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" \
      -o /tmp/fastqc.zip && \
    echo "${FASTQC_SHA256}  /tmp/fastqc.zip" | sha256sum -c - && \
    unzip /tmp/fastqc.zip -d /opt && \
    chmod 755 /opt/FastQC/fastqc && \
    ln -s /opt/FastQC/fastqc /usr/local/bin/fastqc && \
    rm -rf /tmp/fastqc.zip /var/lib/apt/lists/*

ENTRYPOINT ["fastqc"]
