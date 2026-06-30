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
        cmake \
        gfortran \
        libbz2-dev \
        libcurl4-openssl-dev \
        libglpk-dev \
        libhdf5-dev \
        liblzma-dev \
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

# edirect's esearch/efetch are shell scripts that source ecommon.sh from their
# own directory, so they must be run from /root/edirect on PATH -- symlinking the
# individual tools into /usr/local/bin breaks that lookup.
RUN sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
ENV PATH="/root/edirect:${PATH}"

WORKDIR /opt/microsuite
COPY . /opt/microsuite

RUN python3 -m venv /opt/venv \
    && pip install --upgrade pip \
    && pip install ".[biom,qza]"

RUN Rscript -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org"); repos <- c(adw96 = "https://adw96.r-universe.dev", BiocManager::repositories()); install.packages(c("breakaway", "iNEXT"), repos = repos, dependencies = c("Depends", "Imports", "LinkingTo"), Ncpus = max(1, parallel::detectCores() - 1)); missing <- c("breakaway", "iNEXT")[!vapply(c("breakaway", "iNEXT"), requireNamespace, logical(1), quietly = TRUE)]; if (length(missing)) stop("Missing R packages after install: ", paste(missing, collapse = ", "))'

COPY cloud/prjna321534/run_prjna321534.sh /usr/local/bin/run_prjna321534
COPY cloud/prjna321534/run_prjna321534_microsuite.sh /usr/local/bin/run_prjna321534_microsuite
COPY cloud/prjna321534/alpha_breakaway_inext.R /usr/local/bin/alpha_breakaway_inext.R
RUN chmod +x /usr/local/bin/run_prjna321534 /usr/local/bin/run_prjna321534_microsuite

# Default entrypoint runs the raw-tools pipeline; override the command (e.g.
# `docker run IMAGE run_prjna321534_microsuite`) to drive the microsuite CLI.
ENTRYPOINT ["run_prjna321534"]
