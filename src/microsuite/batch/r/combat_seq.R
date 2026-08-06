#!/usr/bin/env Rscript
# sva::ComBat_seq: negative-binomial batch adjustment returning integer counts.
# Usage: combat_seq.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: combat_seq.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(sva)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

meta <- meta[colnames(counts), , drop = FALSE]

# ComBat_seq models counts, so a non-integer matrix means the caller fed it
# something already normalized. Fail rather than round silently.
if (any(abs(counts - round(counts)) > 1e-8)) {
  stop("ComBat_seq requires integer counts; this table holds non-integer values.")
}
storage.mode(counts) <- "integer"

batch <- factor(meta[[params$batch]])

covar_mod <- NULL
if (length(params$covariates) > 0) {
  covariates <- as.character(params$covariates)
  design <- meta[, covariates, drop = FALSE]
  for (name in covariates) {
    if (is.character(design[[name]])) design[[name]] <- factor(design[[name]])
  }
  covar_mod <- model.matrix(~., data = design)
}

adjusted <- ComBat_seq(counts = counts, batch = batch, group = NULL, covar_mod = covar_mod)

out <- data.frame(feature_id = rownames(adjusted), adjusted, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
