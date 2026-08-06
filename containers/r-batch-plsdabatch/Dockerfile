FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-plsdabatch"
LABEL org.opencontainers.image.description="PLSDA-batch CLR-space, supervised batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: PLSDAbatch (GitHub, pinned by commit) (+ mixOmics, jsonlite)
#
# PLSDAbatch's DESCRIPTION at the pinned commit declares Depends: R (>= 4.5.0)
# and Imports: ggplot2, ggpubr, grid, gridExtra, lmerTest, mixOmics,
# performance, scales, stats, Rdpack. grid/stats ship with R itself; the rest
# resolve on conda-forge/bioconda under the usual r-* names, and bioconda
# carries an R-4.5-built bioconductor-mixomics, so r-base is pinned to the
# 4.5 series (unlike the other R batch-correction images here, which use
# 4.3) to satisfy PLSDAbatch's own version floor.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.5,<4.6" \
        r-jsonlite \
        r-remotes \
        bioconductor-mixomics \
        r-ggplot2 \
        r-ggpubr \
        r-gridextra \
        r-lmertest \
        r-performance \
        r-scales \
        r-rdpack && \
    micromamba clean -a -y

# Pinned to a commit, never a branch: reproducibility across rebuilds.
ARG PLSDABATCH_SHA=24f8f48ba77c81c156b0633f4593127bcaf9408e
RUN Rscript -e "remotes::install_github('EvaYiwenWang/PLSDAbatch', ref = '${PLSDABATCH_SHA}', upgrade = 'never')" && \
    Rscript -e "stopifnot(requireNamespace('PLSDAbatch', quietly = TRUE))"

COPY src/microsuite/batch/r/plsda_batch.R /opt/microsuite/plsda_batch.R
COPY containers/r-batch-plsdabatch/smoke/ /opt/microsuite/smoke/

# Build-time smoke: a real correction. This is a supervised backend, so the
# smoke params carry a target. Output is CLR log-ratios (signed,
# non-integer), so unlike combat_seq.R/conqur.R this smoke does NOT assert
# integrality of the output -- that assertion would be wrong for this
# backend and would fail the build.
RUN set -eux; \
    Rscript /opt/microsuite/plsda_batch.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
