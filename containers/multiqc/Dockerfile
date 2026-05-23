FROM python:3.11-slim

LABEL org.opencontainers.image.title="multiqc"
LABEL org.opencontainers.image.description="MultiQC aggregate quality-report backend"

# Expected commands: multiqc
RUN python -m pip install --no-cache-dir multiqc

ENTRYPOINT ["multiqc"]
