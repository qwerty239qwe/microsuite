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
features <- as.data.frame(t(counts), check.names = FALSE)

invisible(maaslin3::maaslin3(
  input_data = features,
  input_metadata = metadata,
  output = output_dir,
  formula = stats::as.formula(params$formula),
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
