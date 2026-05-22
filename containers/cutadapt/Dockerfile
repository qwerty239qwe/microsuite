FROM python:3.12-slim

ARG CUTADAPT_VERSION=4.9

LABEL org.opencontainers.image.title="cutadapt"
LABEL org.opencontainers.image.description="Cutadapt adapter and primer trimming backend"

# Expected commands: cutadapt
RUN pip install --no-cache-dir "cutadapt==${CUTADAPT_VERSION}"

ENTRYPOINT ["cutadapt"]
