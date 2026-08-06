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

# Pinned to a commit, never a branch: `git ls-remote
# https://github.com/BoYuan07/MetaDICT HEAD` resolves the current tip of the
# default branch to this SHA at the time this Dockerfile was written -- that
# confirms what the tip pointed to then, not that this object is otherwise
# distinguished (e.g. a tagged release); it will keep resolving to this same
# commit going forward regardless of what the branch tip does next, which is
# the property pinning by commit (rather than by branch name) is for. The
# repository path itself (BoYuan07/MetaDICT, not the plan's original
# wangyf1996/MetaDICT, which 404s) was corrected in Task 4 and verified twice
# against the repo README before this task began.
ARG METADICT_SHA=5b052877328c05e7337e4ce2789a9f48fdecbd9b
RUN Rscript -e "remotes::install_github('BoYuan07/MetaDICT', ref = '${METADICT_SHA}', upgrade = 'never')" && \
    Rscript -e "stopifnot(requireNamespace('MetaDICT', quietly = TRUE))"

COPY src/microsuite/batch/r/metadict.R /opt/microsuite/metadict.R
COPY containers/r-batch-metadict/smoke/ /opt/microsuite/smoke/

# Build-time smoke: a real correction. MetaDICT declares `relative` as its
# output type (see backends.py); unlike combat_seq.R/conqur.R this is not an
# integer-output assertion but a proportions assertion -- metadict.R
# TSS-normalizes MetaDICT's count-like `$count` output before writing, so
# every sample column must sum to ~1. This run also exercises metadict.R's
# inferred, undocumented fallback (a uniform taxon distance matrix), since
# the smoke params.json supplies no distance_matrix_path -- see the header
# comment in metadict.R.
RUN set -eux; \
    Rscript /opt/microsuite/metadict.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    Rscript -e "x <- read.delim('/tmp/microsuite-smoke-out.tsv', row.names=1); \
                stopifnot(all(abs(colSums(as.matrix(x)) - 1) < 1e-6))"; \
    rm -f /tmp/microsuite-smoke-out.tsv

# Second build-time smoke: metadict.R sorts samples into contiguous batch
# blocks before calling MetaDICT and restores the caller's original column
# order afterward, to neutralize an undocumented column-mislabeling
# behaviour found by reading R/MetaDICT.R (see metadict.R's header
# comment). smoke/counts.tsv is already batch-contiguous (s0-s5=A,
# s6-s11=B), so that first smoke run alone exercises the sort/restore as a
# no-op and would not catch a broken fix. smoke/counts_interleaved.tsv holds
# the identical values, permuted so batches interleave column-by-column
# (s0,s6,s1,s7,...); metadata.tsv is reused unchanged, since metadict.R
# re-aligns metadata to the count columns by name regardless of the
# metadata file's own row order. If the sort/restore is correct, both runs
# feed MetaDICT the same batch-sorted matrix internally and must therefore
# agree, sample-for-sample, once both outputs are reordered to a common
# column order -- this is checked exactly (to floating-point tolerance),
# not just by line count.
RUN set -eux; \
    Rscript /opt/microsuite/metadict.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out-contiguous.tsv; \
    Rscript /opt/microsuite/metadict.R \
        /opt/microsuite/smoke/counts_interleaved.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out-interleaved.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out-interleaved.tsv)" -eq 21; \
    Rscript -e " \
        a <- read.delim('/tmp/microsuite-smoke-out-contiguous.tsv', row.names = 1, check.names = FALSE); \
        b <- read.delim('/tmp/microsuite-smoke-out-interleaved.tsv', row.names = 1, check.names = FALSE); \
        stopifnot(setequal(colnames(a), colnames(b))); \
        stopifnot(setequal(rownames(a), rownames(b))); \
        b <- b[rownames(a), colnames(a)]; \
        diff <- max(abs(as.matrix(a) - as.matrix(b))); \
        message(sprintf('max abs diff between contiguous and interleaved runs: %g', diff)); \
        stopifnot(diff < 1e-6)"; \
    rm -f /tmp/microsuite-smoke-out-contiguous.tsv /tmp/microsuite-smoke-out-interleaved.tsv

# Third build-time smoke: the zero-covariates path. metadict.R's header
# comment documents that a zero-length `covariates` reaches a crash inside
# MetaDICT itself (rbind of zero-column frames collapsing to zero rows), so
# metadict.R substitutes a constant single-level shim column as the sole
# covariate whenever params.covariates is empty -- an inference, not
# something the docs recommend. This is that branch's only exercise, real
# or otherwise, so it runs here rather than waiting for a user to hit it.
RUN set -eux; \
    Rscript /opt/microsuite/metadict.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params_no_covariates.json \
        /tmp/microsuite-smoke-out-no-covariates.tsv; \
    test -s /tmp/microsuite-smoke-out-no-covariates.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out-no-covariates.tsv)" -eq 21; \
    rm -f /tmp/microsuite-smoke-out-no-covariates.tsv

ENTRYPOINT ["Rscript"]
