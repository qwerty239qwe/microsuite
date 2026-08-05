FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-ecology"
LABEL org.opencontainers.image.description="vegan adonis2 and anosim2 beta-diversity backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: vegan
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-vegan && \
    micromamba clean -a -y

COPY src/microsuite/diversity/r/beta_significance.R /opt/microsuite/beta_significance.R
COPY containers/r-ecology/smoke/ /opt/microsuite/smoke/

# Build-time smoke: execute both supported analysis paths on a small real
# distance matrix so a successful package install cannot hide a broken image.
RUN set -eux; \
    Rscript /opt/microsuite/beta_significance.R \
        /opt/microsuite/smoke/distance.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        adonis2 group '' 19 0 /tmp/microsuite-adonis2.tsv; \
    test -s /tmp/microsuite-adonis2.tsv; \
    Rscript /opt/microsuite/beta_significance.R \
        /opt/microsuite/smoke/distance.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        anosim2 group '' 19 0 /tmp/microsuite-anosim2.tsv; \
    test -s /tmp/microsuite-anosim2.tsv; \
    rm -f /tmp/microsuite-adonis2.tsv /tmp/microsuite-anosim2.tsv

ENTRYPOINT ["Rscript"]
