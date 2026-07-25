FROM python:3.11-slim

LABEL org.opencontainers.image.title="microsuite"
LABEL org.opencontainers.image.description="Python CLI and SDK layer for the microsuite toolbox"

# Expected commands: microsuite, uv
# libgomp1 supports downloaded BLAST+ binaries used by FOMC workflows.
RUN apt-get update && \
    apt-get install -y --no-install-recommends ca-certificates libgomp1 procps && \
    rm -rf /var/lib/apt/lists/* && \
    python -m pip install --no-cache-dir uv

WORKDIR /opt/microsuite
COPY . /opt/microsuite
# The command is also used by workflow engines that set an arbitrary host UID.
# Packages are immutable at runtime, but must be readable and executable.
RUN uv sync && chmod -R a+rX /opt/microsuite

ENTRYPOINT ["uv", "run", "microsuite"]
