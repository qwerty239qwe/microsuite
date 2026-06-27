#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: breakaway_alpha.R counts.tsv output.tsv", call. = FALSE)
}

counts_path <- args[[1]]
output_path <- args[[2]]

if (!requireNamespace("breakaway", quietly = TRUE)) {
  stop("R package 'breakaway' is required.", call. = FALSE)
}

counts <- read.delim(counts_path, check.names = FALSE, row.names = 1)

extract_estimate <- function(fit) {
  fields <- unclass(fit)
  estimate <- NA_real_
  error <- NA_real_
  for (name in c("estimate", "estimated_richness", "richness")) {
    if (!is.null(fields[[name]])) {
      estimate <- as.numeric(fields[[name]])[[1]]
      break
    }
  }
  for (name in c("error", "std_error", "se", "standard_error")) {
    if (!is.null(fields[[name]])) {
      error <- as.numeric(fields[[name]])[[1]]
      break
    }
  }
  data.frame(estimate = estimate, error = error, stringsAsFactors = FALSE)
}

rows <- list()
for (sample_id in colnames(counts)) {
  values <- counts[[sample_id]]
  values <- values[values > 0]
  result <- tryCatch(
    {
      fit <- breakaway::breakaway(values)
      extracted <- extract_estimate(fit)
      data.frame(
        sample_id = sample_id,
        method = "breakaway",
        estimate = extracted$estimate,
        error = extracted$error,
        status = "ok",
        message = "",
        stringsAsFactors = FALSE
      )
    },
    error = function(err) {
      data.frame(
        sample_id = sample_id,
        method = "breakaway",
        estimate = NA_real_,
        error = NA_real_,
        status = "error",
        message = conditionMessage(err),
        stringsAsFactors = FALSE
      )
    }
  )
  rows[[length(rows) + 1]] <- result
}

write.table(
  do.call(rbind, rows),
  file = output_path,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
