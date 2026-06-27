#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(
    paste(
      "Usage: inext_alpha.R counts.tsv output.tsv q_csv datatype knots",
      "se conf nboot"
    ),
    call. = FALSE
  )
}

counts_path <- args[[1]]
output_path <- args[[2]]
q <- as.numeric(strsplit(args[[3]], ",", fixed = TRUE)[[1]])
datatype <- args[[4]]
knots <- as.integer(args[[5]])
se <- tolower(args[[6]]) == "true"
conf <- as.numeric(args[[7]])
nboot <- as.integer(args[[8]])

if (!requireNamespace("iNEXT", quietly = TRUE)) {
  stop("R package 'iNEXT' is required.", call. = FALSE)
}

counts <- read.delim(counts_path, check.names = FALSE, row.names = 1)
sample_list <- lapply(counts, function(values) values[values > 0])

rows <- list()
for (sample_id in names(sample_list)) {
  result <- tryCatch(
    {
      fit <- iNEXT::iNEXT(
        list(sample = sample_list[[sample_id]]),
        q = q,
        datatype = datatype,
        knots = knots,
        se = se,
        conf = conf,
        nboot = nboot
      )
      asy <- fit$AsyEst
      asy$sample_id <- sample_id
      asy$method <- "inext"
      asy$status <- "ok"
      asy$message <- ""
      asy
    },
    error = function(err) {
      data.frame(
        sample_id = sample_id,
        method = "inext",
        Observed = NA_real_,
        Estimator = NA_real_,
        s.e. = NA_real_,
        LCL = NA_real_,
        UCL = NA_real_,
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
