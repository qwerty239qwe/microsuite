args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: maaslin3.R counts.tsv metadata.tsv params.json output_dir")
}

counts_path <- args[[1]]
metadata_path <- args[[2]]
params_path <- args[[3]]
output_dir <- args[[4]]

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("R package 'jsonlite' is required.")
}
if (!requireNamespace("maaslin3", quietly = TRUE)) {
  stop("R package 'maaslin3' is required. Install it with BiocManager::install('maaslin3').")
}

params <- jsonlite::fromJSON(params_path, simplifyVector = TRUE)
required <- c("formula", "normalization", "transform", "min_prevalence", "min_abundance")
missing_params <- setdiff(required, names(params))
if (length(missing_params) > 0) {
  stop(paste("Missing MaAsLin 3 parameters:", paste(missing_params, collapse = ", ")))
}

counts <- read.delim(counts_path, row.names = 1, check.names = FALSE)
metadata <- read.delim(metadata_path, row.names = 1, check.names = FALSE)
missing_metadata <- setdiff(colnames(counts), rownames(metadata))
if (length(missing_metadata) > 0) {
  stop(paste("Samples missing from metadata:", paste(missing_metadata, collapse = ", ")))
}
metadata <- metadata[colnames(counts), , drop = FALSE]

formula <- stats::as.formula(params$formula)
term_labels <- attr(stats::terms(formula), "term.labels")
special_effect <- grepl("^(group|ordered|strata)\\s*\\(", term_labels)
fixed_labels <- term_labels[!special_effect & !grepl("\\|", term_labels)]
fixed_formula <- if (length(fixed_labels) == 0L) {
  stats::as.formula("~ 1")
} else {
  stats::as.formula(paste("~", paste(fixed_labels, collapse = " + ")))
}
fixed_variables <- all.vars(fixed_formula)

# TSV input cannot retain pandas categorical dtypes. MaAsLin 3 consequently
# receives multi-level categorical fixed effects as character vectors and asks
# for a reference. Promote only modeled fixed-effect character columns to
# factors. Without an explicit reference the first sorted level is the baseline.
reference_pairs <- list()
if (!is.null(params$reference) && length(params$reference) > 0) {
  reference <- as.character(params$reference)
  if (length(reference) != 1L || !nzchar(trimws(reference))) {
    stop("MaAsLin 3 reference must be one non-empty string")
  }
  for (pair in strsplit(reference, ";", fixed = TRUE)[[1]]) {
    parts <- strsplit(trimws(pair), ",", fixed = TRUE)[[1]]
    if (length(parts) != 2L || any(!nzchar(trimws(parts)))) {
      stop("Invalid MaAsLin 3 reference; expected 'column,level' pairs separated by ';'")
    }
    column <- trimws(parts[[1]])
    level <- trimws(parts[[2]])
    if (!is.null(reference_pairs[[column]])) {
      stop(paste0("MaAsLin 3 reference column is specified more than once: ", column))
    }
    if (!(column %in% colnames(metadata))) {
      stop(paste0("MaAsLin 3 reference column not found in sample metadata: ", column))
    }
    if (!(column %in% fixed_variables)) {
      stop(paste0("MaAsLin 3 reference column is not a fixed effect in the formula: ", column))
    }
    values <- metadata[[column]]
    if (!(is.character(values) || is.factor(values) || is.logical(values))) {
      stop(paste0("MaAsLin 3 reference column must be categorical: ", column))
    }
    if (!(level %in% unique(as.character(values)))) {
      stop(paste0("Reference level '", level, "' not found in metadata column '", column, "'"))
    }
    reference_pairs[[column]] <- level
  }
}

for (column in intersect(fixed_variables, colnames(metadata))) {
  values <- metadata[[column]]
  wanted <- reference_pairs[[column]]
  if (is.character(values) || is.logical(values)) {
    levels_sorted <- sort(unique(as.character(values)))
    if (!is.null(wanted)) levels_sorted <- c(wanted, setdiff(levels_sorted, wanted))
    metadata[[column]] <- factor(as.character(values), levels = levels_sorted)
  }
  if (is.factor(metadata[[column]]) && nlevels(metadata[[column]]) > 2L) {
    message(paste0(
      "microsuite: '", column, "' reference level = '", levels(metadata[[column]])[[1]], "'"
    ))
  }
}

features <- as.data.frame(t(counts), check.names = FALSE)

invisible(maaslin3::maaslin3(
  input_data = features,
  input_metadata = metadata,
  output = output_dir,
  formula = formula,
  reference = if (is.null(params$reference)) NULL else as.character(params$reference),
  normalization = params$normalization,
  transform = params$transform,
  min_prevalence = params$min_prevalence,
  min_abundance = params$min_abundance,
  plot_summary_plot = FALSE,
  plot_associations = FALSE,
  save_models = FALSE,
  cores = 1,
  verbosity = "WARN"
))

all_results_path <- file.path(output_dir, "all_results.tsv")
if (!file.exists(all_results_path)) {
  stop("MaAsLin 3 did not produce all_results.tsv")
}
results <- read.delim(all_results_path, check.names = FALSE)
if (!("model" %in% colnames(results))) {
  stop("MaAsLin 3 all_results.tsv is missing the model column")
}

model <- tolower(as.character(results$model))
abundance <- results[model == "abundance", , drop = FALSE]
prevalence <- results[model == "prevalence", , drop = FALSE]
write.table(
  abundance,
  file = file.path(output_dir, "abundance_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  prevalence,
  file = file.path(output_dir, "prevalence_results.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
invisible(file.copy(params_path, file.path(output_dir, "microsuite_params.json"), overwrite = TRUE))
