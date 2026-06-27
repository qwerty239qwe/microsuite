#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: alpha_breakaway_inext.R otu_table.tsv output_dir", call. = FALSE)
}

otu_path <- args[[1]]
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("breakaway", quietly = TRUE)) {
  stop("R package 'breakaway' is required.", call. = FALSE)
}
if (!requireNamespace("iNEXT", quietly = TRUE)) {
  stop("R package 'iNEXT' is required.", call. = FALSE)
}

otu <- read.delim(otu_path, check.names = FALSE, row.names = 1, comment.char = "")

extract_breakaway <- function(fit) {
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

breakaway_rows <- list()
inext_rows <- list()

sample_list <- lapply(otu, function(values) as.numeric(values[values > 0]))

empty_inext_row <- function(sample_id, status, message) {
  data.frame(
    sample_id = sample_id,
    method = "inext",
    Assemblage = NA_character_,
    Diversity = NA_character_,
    Observed = NA_real_,
    Estimator = NA_real_,
    s.e. = NA_real_,
    LCL = NA_real_,
    UCL = NA_real_,
    status = status,
    message = message,
    stringsAsFactors = FALSE
  )
}

normalize_inext <- function(asy, sample_id) {
  required <- empty_inext_row(sample_id, "ok", "")
  asy$sample_id <- sample_id
  asy$method <- "inext"
  asy$status <- "ok"
  asy$message <- ""
  for (name in names(required)) {
    if (!name %in% names(asy)) {
      asy[[name]] <- required[[name]]
    }
  }
  asy[names(required)]
}

for (sample_id in names(sample_list)) {
  values <- sample_list[[sample_id]]

  breakaway_rows[[length(breakaway_rows) + 1]] <- tryCatch(
    {
      fit <- breakaway::breakaway(values)
      extracted <- extract_breakaway(fit)
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

  inext_rows[[length(inext_rows) + 1]] <- tryCatch(
    {
      fit <- iNEXT::iNEXT(
        list(sample = values),
        q = c(0, 1, 2),
        datatype = "abundance",
        knots = 40,
        se = TRUE,
        conf = 0.95,
        nboot = 50
      )
      normalize_inext(fit$AsyEst, sample_id)
    },
    error = function(err) {
      empty_inext_row(sample_id, "error", conditionMessage(err))
    }
  )
}

write.table(
  do.call(rbind, breakaway_rows),
  file = file.path(output_dir, "breakaway_alpha.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  do.call(rbind, inext_rows),
  file = file.path(output_dir, "inext_alpha.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
