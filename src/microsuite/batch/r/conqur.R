#!/usr/bin/env Rscript
# ConQuR: conditional quantile regression batch removal, returning counts.
#
# UNVERIFIED AT RUNTIME: this signature is taken from the published man page
# (man/ConQuR.Rd at the pinned commit, https://github.com/wdl2459/ConQuR) --
# there is no container engine available in the environment this script was
# written in, so it has never been executed against the real ConQuR package.
# The container's build-time smoke test (see containers/r-batch-conqur/Dockerfile)
# is this script's first real execution. If that smoke fails, re-check this
# script against `Rscript -e "args(ConQuR::ConQuR)"` inside the built image.
#
# Usage: conqur.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: conqur.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(ConQuR)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

# ConQuR::ConQuR takes tax_tab as samples (row) by taxa (col); microsuite
# writes features as rows, samples as columns, so transpose in and back out.
tax_tab <- as.data.frame(t(counts))
batchid <- factor(meta[[params$batch]])

# ConQuR requires a reference batch and documents no default. Default to the
# first level here so runs are reproducible rather than dependent on factor
# ordering elsewhere.
batch_ref <- if (!is.null(params$batch_ref)) as.character(params$batch_ref) else levels(batchid)[1]

if (length(params$covariates) > 0) {
  covariates <- meta[, as.character(params$covariates), drop = FALSE]
  for (name in names(covariates)) {
    if (is.character(covariates[[name]])) covariates[[name]] <- factor(covariates[[name]])
  }
} else {
  # ConQuR requires a covariate data.frame; an intercept-only frame is the
  # no-covariate case.
  covariates <- data.frame(intercept_only = factor(rep("a", nrow(tax_tab))))
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
