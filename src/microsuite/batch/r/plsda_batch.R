#!/usr/bin/env Rscript
# PLSDA-batch: subtracts batch-associated latent components in CLR space.
# Output is CLR log-ratios, NOT counts or relative abundances.
#
# UNVERIFIED AT RUNTIME: this signature is taken from the published man page
# (man/PLSDA_batch.Rd at the pinned commit,
# https://github.com/EvaYiwenWang/PLSDAbatch) -- there is no container engine
# available in the environment this script was written in, so it has never
# been executed against the real PLSDAbatch package. The container's
# build-time smoke test (see containers/r-batch-plsdabatch/Dockerfile) is
# this script's first real execution. If that smoke fails, re-check this
# script against `Rscript -e "args(PLSDAbatch::PLSDA_batch)"` inside the
# built image.
#
# man/PLSDA_batch.Rd confirms:
#   PLSDA_batch(X, Y.trt = NULL, Y.bat, ncomp.trt = 2, ncomp.bat = 2, ...)
#   X: "a numeric matrix as an explanatory matrix" -- the man page's own
#   example builds X from a CLR-transformed assay (`Clr_value`) indexed the
#   same way as the per-sample Y.trt/Y.bat vectors (via rownames() of the
#   source object), i.e. X is samples-as-rows, CLR-transformed already, not
#   raw counts. Return value: a list including `X.nobatch`, "the batch
#   corrected matrix with the same dimension as the input matrix" -- matches
#   the plan's assumption exactly.
#
# Usage: plsda_batch.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: plsda_batch.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(mixOmics)
  library(PLSDAbatch)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

if (is.null(params$target) || is.na(params$target)) {
  stop("plsda-batch is supervised and requires a target column.")
}

# PLSDA_batch's X is samples-as-rows; microsuite writes features as rows,
# samples as columns, so transpose in and back out (as conqur.R does).
#
# man/PLSDA_batch.Rd does not itself specify how to get from raw counts to
# the CLR input it expects; it only shows a pre-computed `Clr_value` assay
# in its example. The CLR transform used here is mixOmics::logratio.transfo,
# documented in mixOmics man/logratio-transformations.Rd as
# `logratio.transfo(X, logratio = c("none", "CLR", "ILR"), offset = 0)`,
# where X is samples-as-rows and `offset` "is used to shift the values away
# from 0, as commonly done with counts data" for CLR/ILR. That confirms the
# offset mechanism and its purpose, but not a value to use for count data:
# the parameter's own default is 0 (a no-op), and the man page's CLR example
# applies no offset at all (its input is already non-zero after TSS).
# INFERRED, NOT DOCUMENTED: an offset of 1 is used here, the common
# convention for zero-laden count data (matching combat_seq.R/conqur.R's
# treatment of zeros elsewhere in this codebase), not something
# logratio-transformations.Rd itself recommends.
X <- t(counts)
X <- logratio.transfo(X = X, logratio = "CLR", offset = 1)
class(X) <- "matrix"

Y.trt <- factor(meta[[params$target]])
Y.bat <- factor(meta[[params$batch]])

ncomp_trt <- if (!is.null(params$ncomp_trt)) as.integer(params$ncomp_trt) else 2L
ncomp_bat <- if (!is.null(params$ncomp_bat)) as.integer(params$ncomp_bat) else 2L

fit <- PLSDA_batch(
  X = X,
  Y.trt = Y.trt,
  Y.bat = Y.bat,
  ncomp.trt = ncomp_trt,
  ncomp.bat = ncomp_bat
)
adjusted <- fit$X.nobatch

# Back to features-as-rows for the caller. Output is CLR log-ratios (signed,
# non-integer); no integrality check applies here, unlike the counts-typed
# backends.
out_matrix <- t(as.matrix(adjusted))
out <- data.frame(feature_id = rownames(out_matrix), out_matrix, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
