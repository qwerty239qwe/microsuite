FROM ubuntu:24.04

LABEL org.opencontainers.image.title="microsuite-prjna321534-alpha"
LABEL org.opencontainers.image.description="SRA/VSEARCH/Breakaway/iNEXT runner for BioProject PRJNA321534 alpha-diversity smoke and production runs"

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/opt/venv/bin:${PATH}"

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        default-jre-headless \
        git \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        python3 \
        python3-pip \
        python3-venv \
        r-base \
        r-base-dev \
        sra-toolkit \
        vsearch \
    && rm -rf /var/lib/apt/lists/*

RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)" \
    && ln -s /root/edirect/esearch /usr/local/bin/esearch \
    && ln -s /root/edirect/efetch /usr/local/bin/efetch \
    && ln -s /root/edirect/xtract /usr/local/bin/xtract

WORKDIR /opt/microsuite
COPY . /opt/microsuite

RUN python3 -m venv /opt/venv \
    && pip install --upgrade pip \
    && pip install ".[biom,qza]"

RUN Rscript -e 'install.packages(c("breakaway", "iNEXT"), repos = "https://cloud.r-project.org")'

COPY cloud/prjna321534/run_prjna321534.sh /usr/local/bin/run_prjna321534
COPY cloud/prjna321534/alpha_breakaway_inext.R /usr/local/bin/alpha_breakaway_inext.R

ENTRYPOINT ["run_prjna321534"]
