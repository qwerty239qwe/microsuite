args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop("usage: tax4fun2.R <seqs.fasta> <otu-table.tsv> <reference-dir> <output-dir> <threads> <database-mode> <min-identity> <normalize-pathways>", call. = FALSE)
}

seqs <- args[[1]]
otu_table <- args[[2]]
reference_dir <- args[[3]]
output_dir <- args[[4]]
threads <- as.integer(args[[5]])
database_mode <- args[[6]]
min_identity <- as.numeric(args[[7]])
normalize_pathways <- tolower(args[[8]]) %in% c("1", "true", "t", "yes")

if (!requireNamespace("Tax4Fun2", quietly = TRUE)) {
  stop("Tax4Fun2 R package is required.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
temp_dir <- file.path(output_dir, paste0("tax4fun2_", database_mode))

Tax4Fun2::runRefBlast(
  path_to_otus = seqs,
  path_to_reference_data = reference_dir,
  path_to_temp_folder = temp_dir,
  database_mode = database_mode,
  use_force = TRUE,
  num_threads = threads
)

Tax4Fun2::makeFunctionalPrediction(
  path_to_otu_table = otu_table,
  path_to_reference_data = reference_dir,
  path_to_temp_folder = temp_dir,
  database_mode = database_mode,
  normalize_by_copy_number = TRUE,
  min_identity_to_reference = min_identity,
  normalize_pathways = normalize_pathways
)
