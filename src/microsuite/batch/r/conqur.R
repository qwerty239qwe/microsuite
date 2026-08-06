#!/usr/bin/env Rscript
# ConQuR: conditional quantile regression batch removal, returning counts.
#
# Signature source: the published man page (man/ConQuR.Rd at the pinned commit,
# https://github.com/wdl2459/ConQuR). It was written without a container engine
# available, so the r-batch-conqur build-time smoke was its first execution.
# That smoke has since run and caught a missing `foreach` attachment (see the
# library() block below). Whether the rest of the call is right is decided by
# the next heavy-image build, not by this comment.
#
# Usage: conqur.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: conqur.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  # ConQuR calls `foreach(...) %do% {...}` but its NAMESPACE only has
  # import(doParallel) -- it never imports foreach's `%do%` operator, and
  # importing doParallel does not re-export it. Inside ConQuR's namespace the
  # lookup for `%do%` therefore falls through to the search path, so foreach
  # must be ATTACHED here or every call fails with:
  #   Error in foreach(...) %do% { : could not find function "%do%"
  # Verified against DESCRIPTION and NAMESPACE at the pinned commit
  # c7a88794efd4ecfe4d96988dceeec3b410222e48, and reproduced by the
  # r-batch-conqur build-time smoke before this line was added.
  library(foreach)
  library(doParallel)
  library(ConQuR)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

# Unlike ComBat_seq (a negative-binomial count model), ConQuR fits
# conditional quantile regressions over the observed values; nothing in
# man/ConQuR.Rd requires integer input, so this script does not reject
# non-integer counts the way combat_seq.R does. The build-time smoke still
# checks the *output* is integer, since that is the capability microsuite
# declares for this backend.
#
# ConQuR::ConQuR takes tax_tab as samples (row) by taxa (col); microsuite
# writes features as rows, samples as columns, so transpose in and back out.
tax_tab <- as.data.frame(t(counts))
batchid <- factor(meta[[params$batch]])

# ConQuR requires a reference batch and documents no default. Default to the
# first level here so runs are reproducible rather than dependent on factor
# ordering elsewhere.
batch_ref <- if (!is.null(params$batch_ref)) as.character(params$batch_ref) else levels(batchid)[1]

if (length(params$covariates) == 0) {
  # Discovered by execution, not from the docs: ConQuR builds a model matrix
  # over `covariates`, so an empty set leaves nothing to condition on. The
  # earlier intercept-only shim (a single-level factor) reached
  # model.matrix and died with "contrasts can be applied only to factors with
  # 2 or more levels". ConQuR is conditional by construction, so refusing is
  # correct: the caller must name the variable the correction has to preserve.
  stop(
    "ConQuR requires at least one covariate: it removes batch effects while ",
    "holding the named variables fixed, and its design matrix is degenerate ",
    "without one. Pass --covariates naming the biological variable of interest."
  )
}
covariates <- meta[, as.character(params$covariates), drop = FALSE]
for (name in names(covariates)) {
  if (is.character(covariates[[name]])) covariates[[name]] <- factor(covariates[[name]])
}

adjusted <- ConQuR(
  tax_tab = tax_tab,
  batchid = batchid,
  covariates = covariates,
  batch_ref = batch_ref
)

# ConQuR returns samples (row) by taxa (col); back to features-as-rows for
# the caller.
out_matrix <- t(as.matrix(adjusted))
out <- data.frame(feature_id = rownames(out_matrix), out_matrix, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
