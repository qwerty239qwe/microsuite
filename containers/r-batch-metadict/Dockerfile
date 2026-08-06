FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-metadict"
LABEL org.opencontainers.image.description="MetaDICT covariate-balancing plus shared-dictionary-learning batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: MetaDICT (GitHub, pinned by commit) (+ jsonlite)
#
# MetaDICT's Imports (DESCRIPTION at the pinned commit): stats, RANN,
# igraph, vegan, edgeR, ecodist, ggplot2, viridis, ggpubr, ape, cluster,
# matrixStats. Every one of those resolves on conda-forge/bioconda under the
# usual r-*/bioconductor-* names; unlike ConQuR's cqrReg, nothing here is
# CRAN-only or requires a compiler toolchain, so there is no separate
# install-from-source step.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        r-remotes \
        r-rann \
        r-igraph \
        r-vegan \
        bioconductor-edger \
        r-ecodist \
        r-ggplot2 \
        r-viridis \
        r-ggpubr \
        r-ape \
        r-cluster \
        r-matrixstats && \
    micromamba clean -a -y

# Pinned to a commit, never a branch: this SHA was confirmed via
# `git ls-remote https://github.com/BoYuan07/MetaDICT HEAD`. The repository
# path itself (BoYuan07/MetaDICT, not the plan's original wangyf1996/MetaDICT,
# which 404s) was corrected in Task 4 and verified twice against the repo
# README before this task began.
ARG METADICT_SHA=5b052877328c05e7337e4ce2789a9f48fdecbd9b
RUN Rscript -e "remotes::install_github('BoYuan07/MetaDICT', ref = '${METADICT_SHA}', upgrade = 'never')" && \
    Rscript -e "stopifnot(requireNamespace('MetaDICT', quietly = TRUE))"

COPY src/microsuite/batch/r/metadict.R /opt/microsuite/metadict.R
COPY containers/r-batch-metadict/smoke/ /opt/microsuite/smoke/

# Build-time smoke: a real correction. MetaDICT declares `relative` as its
# output type (see backends.py), not counts, so unlike combat_seq.R/conqur.R
# there is no integer-output assertion here -- only the line-count check.
# This run also exercises metadict.R's inferred, undocumented fallback (a
# uniform taxon distance matrix), since the smoke params.json supplies no
# distance_matrix_path -- see the header comment in metadict.R.
RUN set -eux; \
    Rscript /opt/microsuite/metadict.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    rm -f /tmp/microsuite-smoke-out.tsv

ENTRYPOINT ["Rscript"]
