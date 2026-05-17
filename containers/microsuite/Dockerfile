FROM python:3.11-slim

LABEL org.opencontainers.image.title="microsuite"
LABEL org.opencontainers.image.description="Python CLI and SDK layer for the microsuite toolbox"

# Expected commands: microsuite, uv
RUN python -m pip install --no-cache-dir uv

WORKDIR /opt/microsuite
COPY . /opt/microsuite
RUN uv sync

ENTRYPOINT ["uv", "run", "microsuite"]
