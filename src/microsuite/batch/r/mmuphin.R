#!/usr/bin/env Rscript
# MMUPHin adjust_batch: ComBat extended to zero-inflated microbiome profiles.
# Usage: mmuphin.R counts.tsv metadata.tsv params.json corrected.tsv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("usage: mmuphin.R counts.tsv metadata.tsv params.json corrected.tsv")
}
suppressPackageStartupMessages({
  library(jsonlite)
  library(MMUPHin)
})

counts <- as.matrix(read.delim(args[1], row.names = 1, check.names = FALSE))
meta <- read.delim(args[2], row.names = 1, check.names = FALSE)
params <- fromJSON(args[3], simplifyVector = TRUE)

# The caller writes features as rows and samples as columns; align metadata to
# the count columns rather than trusting the two files to share an order.
meta <- meta[colnames(counts), , drop = FALSE]
meta[[params$batch]] <- factor(meta[[params$batch]])

covariates <- if (length(params$covariates) > 0) as.character(params$covariates) else NULL

fit <- adjust_batch(
  feature_abd = counts,
  batch = params$batch,
  covariates = covariates,
  data = meta
)
adjusted <- fit$feature_abd_adj

# adjust_batch()'s documented return is "the same format as feature_abd" --
# i.e. the same scale as its input, which here is raw counts. backends.py
# declares this backend's output as `relative` (proportions summing to ~1
# per sample), so TSS-normalize the corrected matrix here to make that
# label true, rather than passing count-scale numbers downstream under a
# `relative` stamp.
sample_sums <- colSums(adjusted)
adjusted <- sweep(adjusted, 2, sample_sums, "/")

out <- data.frame(feature_id = rownames(adjusted), adjusted, check.names = FALSE)
write.table(out, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
