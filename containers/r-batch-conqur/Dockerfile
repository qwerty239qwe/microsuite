FROM mambaorg/micromamba:1.5.10

LABEL org.opencontainers.image.title="r-batch-conqur"
LABEL org.opencontainers.image.description="ConQuR quantile-regression batch-correction backend for microsuite"

ENV PATH="/opt/conda/bin:${PATH}"

# Expected commands: Rscript
# Expected R package: ConQuR (GitHub, pinned by commit) (+ jsonlite)
#
# ConQuR's Imports (DESCRIPTION at the pinned commit): quantreg, cqrReg,
# glmnet, dplyr, doParallel, gplots, vegan, ade4, compositions, randomForest,
# ROCR, ape, GUniFrac, fastDummies. All of those resolve on conda-forge/
# bioconda under the usual r-* names except cqrReg, which is CRAN-only (not
# packaged for conda-forge or bioconda) and needs compilation (it links
# Rcpp/RcppArmadillo), so it is installed from CRAN in a separate step with a
# conda-forge compiler toolchain.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        "r-base>=4.3,<4.4" \
        r-jsonlite \
        r-remotes \
        r-quantreg \
        r-glmnet \
        r-doparallel \
        r-dplyr \
        r-gplots \
        r-vegan \
        r-ade4 \
        r-compositions \
        r-randomforest \
        r-rocr \
        r-ape \
        r-gunifrac \
        r-fastdummies \
        r-rcpp \
        r-rcpparmadillo \
        compilers && \
    micromamba clean -a -y

# cqrReg is a live CRAN package (not on conda-forge/bioconda) that ConQuR
# imports. It has compiled (Rcpp/RcppArmadillo) code, hence the conda-forge
# "compilers" meta-package above and `micromamba run` here, which activates
# the environment so the compiler wrapper scripts set CC/CXX before R builds
# it from source.
RUN micromamba run -n base Rscript -e "install.packages('cqrReg', repos = 'https://cloud.r-project.org')" && \
    micromamba run -n base Rscript -e "stopifnot(requireNamespace('cqrReg', quietly = TRUE))"

# Pinned to a commit, never a branch: ConQuR has no release tags, so a branch
# reference would let this image change with no change on our side.
ARG CONQUR_SHA=c7a88794efd4ecfe4d96988dceeec3b410222e48
RUN Rscript -e "remotes::install_github('wdl2459/ConQuR', ref = '${CONQUR_SHA}', upgrade = 'never')" && \
    Rscript -e "stopifnot(requireNamespace('ConQuR', quietly = TRUE))"

COPY src/microsuite/batch/r/conqur.R /opt/microsuite/conqur.R
COPY containers/r-batch-conqur/smoke/ /opt/microsuite/smoke/

# Build-time smoke: a real correction, with an integer-output check. ConQuR
# declares `counts` as its output type, so non-integer output would silently
# break every count-requiring downstream method.
RUN set -eux; \
    Rscript /opt/microsuite/conqur.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params.json \
        /tmp/microsuite-smoke-out.tsv; \
    test -s /tmp/microsuite-smoke-out.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out.tsv)" -eq 21; \
    Rscript -e "x <- read.delim('/tmp/microsuite-smoke-out.tsv', row.names=1); \
                stopifnot(all(abs(as.matrix(x) - round(as.matrix(x))) < 1e-8))"; \
    rm -f /tmp/microsuite-smoke-out.tsv

# Second build-time smoke: the zero-covariates fallback in conqur.R (an
# intercept-only, single-level factor) is an inference from undocumented
# behaviour, not something man/ConQuR.Rd specifies -- see the comment in
# conqur.R. This is the branch's only exercise, real or otherwise, so it
# must run here rather than waiting for a user to hit it first.
RUN set -eux; \
    Rscript /opt/microsuite/conqur.R \
        /opt/microsuite/smoke/counts.tsv \
        /opt/microsuite/smoke/metadata.tsv \
        /opt/microsuite/smoke/params_no_covariates.json \
        /tmp/microsuite-smoke-out-no-covariates.tsv; \
    test -s /tmp/microsuite-smoke-out-no-covariates.tsv; \
    test "$(wc -l < /tmp/microsuite-smoke-out-no-covariates.tsv)" -eq 21; \
    Rscript -e "x <- read.delim('/tmp/microsuite-smoke-out-no-covariates.tsv', row.names=1); \
                stopifnot(all(abs(as.matrix(x) - round(as.matrix(x))) < 1e-8))"; \
    rm -f /tmp/microsuite-smoke-out-no-covariates.tsv

ENTRYPOINT ["Rscript"]
