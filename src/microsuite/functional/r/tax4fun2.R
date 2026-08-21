args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
  stop(
    paste(
      "usage: tax4fun2.R <seqs.fasta> <otu-table.tsv> <reference-dir>",
      "<threads> <database-mode> <min-identity> <normalize-pathways> <output-dir>"
    ),
    call. = FALSE
  )
}

# Preserve the original public entrypoint, whose output directory was argument 4.
if (!is.na(suppressWarnings(as.integer(args[[4]])))) {
  seqs <- args[[1]]
  otu_table_path <- args[[2]]
  reference_dir <- args[[3]]
  threads_text <- args[[4]]
  database_mode <- args[[5]]
  min_identity_text <- args[[6]]
  normalize_pathways_text <- args[[7]]
  output_dir <- args[[8]]
} else {
  seqs <- args[[1]]
  otu_table_path <- args[[2]]
  reference_dir <- args[[3]]
  output_dir <- args[[4]]
  threads_text <- args[[5]]
  database_mode <- args[[6]]
  min_identity_text <- args[[7]]
  normalize_pathways_text <- args[[8]]
}

fail <- function(message) stop(message, call. = FALSE)
required_file <- function(path, label) {
  if (!file.exists(path) || dir.exists(path)) fail(paste0(label, " not found: ", path))
}

required_file(seqs, "Representative FASTA")
required_file(otu_table_path, "Tax4Fun2 table")
if (!dir.exists(reference_dir)) fail(paste0("Tax4Fun2 reference directory not found: ", reference_dir))
if (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE, no.. = TRUE)) > 0) {
  fail(paste0("Tax4Fun2 output directory is not empty: ", output_dir))
}

threads <- suppressWarnings(as.integer(threads_text))
if (length(threads) != 1 || is.na(threads) || threads < 1 || as.character(threads) != threads_text) {
  fail("Tax4Fun2 threads must be a positive integer.")
}
if (!database_mode %in% c("Ref99NR", "Ref100NR")) {
  fail("Tax4Fun2 database mode must be Ref99NR or Ref100NR.")
}
min_identity <- suppressWarnings(as.numeric(min_identity_text))
if (length(min_identity) != 1 || !is.finite(min_identity) || min_identity <= 0 || min_identity > 1) {
  fail("Tax4Fun2 min identity must be greater than 0 and less than or equal to 1.")
}
normalized_boolean <- tolower(trimws(normalize_pathways_text))
if (!normalized_boolean %in% c("true", "false", "1", "0", "t", "f", "yes", "no")) {
  fail("Tax4Fun2 normalize-pathways must be a Boolean value.")
}
normalize_pathways <- normalized_boolean %in% c("true", "1", "t", "yes")

if (!requireNamespace("Tax4Fun2", quietly = TRUE)) fail("Tax4Fun2 R package is required.")
if (as.character(utils::packageVersion("Tax4Fun2")) != "1.1.5") {
  fail(paste0(
    "Unsupported Tax4Fun2 version ", utils::packageVersion("Tax4Fun2"),
    "; microsuite requires exactly 1.1.5."
  ))
}
if (!requireNamespace("jsonlite", quietly = TRUE)) fail("jsonlite R package is required.")

blastn <- Sys.which("blastn")
makeblastdb <- Sys.which("makeblastdb")
if (!nzchar(blastn)) fail("BLAST+ blastn was not found on PATH.")
if (!nzchar(makeblastdb)) fail("BLAST+ makeblastdb was not found on PATH.")

ref_fasta <- file.path(reference_dir, database_mode, paste0(database_mode, ".fasta"))
required_files <- c(
  ref_fasta,
  file.path(reference_dir, "KEGG", "ko.txt"),
  file.path(reference_dir, "KEGG", "ko2ptw.txt"),
  file.path(reference_dir, "KEGG", "ptw.txt")
)
for (path in required_files) required_file(path, "Tax4Fun2 reference file")
profile_files <- list.files(
  file.path(reference_dir, database_mode), pattern = "\\.tbl\\.gz$", full.names = TRUE
)
if (length(profile_files) == 0) fail("Tax4Fun2 reference contains no compressed profiles.")

fasta_lines <- readLines(seqs, warn = FALSE)
headers <- grep("^>", fasta_lines, value = TRUE)
if (length(headers) == 0) fail("Representative FASTA contains no records.")
fasta_ids <- sub("[[:space:]].*$", "", substring(headers, 2))
if (any(!nzchar(fasta_ids)) || anyDuplicated(fasta_ids)) {
  fail("Representative FASTA identifiers must be non-empty and unique.")
}

otu_table <- tryCatch(
  read.delim(
    otu_table_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  ),
  error = function(error) fail(paste0("Could not read Tax4Fun2 table: ", conditionMessage(error)))
)
if (ncol(otu_table) < 2 || nrow(otu_table) < 1) {
  fail("Tax4Fun2 table requires feature IDs and at least one non-empty sample column.")
}
feature_ids <- trimws(as.character(otu_table[[1]]))
if (any(!nzchar(feature_ids)) || anyDuplicated(feature_ids)) {
  fail("Tax4Fun2 table feature IDs must be non-empty and unique.")
}
sample_names <- names(otu_table)[-1]
if (any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
  fail("Tax4Fun2 table sample names must be non-empty and unique.")
}
abundance <- suppressWarnings(as.matrix(data.frame(lapply(otu_table[-1], as.numeric))))
if (anyNA(abundance) || any(!is.finite(abundance)) || any(abundance < 0)) {
  fail("Tax4Fun2 abundances must be numeric, finite, and non-negative.")
}
if (any(colSums(abundance) <= 0)) fail("Tax4Fun2 samples must have positive total abundance.")
if (!setequal(feature_ids, fasta_ids)) {
  fail("Tax4Fun2 table and representative FASTA feature IDs must match exactly.")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
upstream_dir <- file.path(output_dir, "upstream")
dir.create(upstream_dir, recursive = TRUE, showWarnings = FALSE)
db_prefix <- file.path(upstream_dir, "reference")
makeblastdb_stdout <- file.path(upstream_dir, "makeblastdb.stdout.log")
makeblastdb_stderr <- file.path(upstream_dir, "makeblastdb.stderr.log")
status <- system2(
  makeblastdb,
  args = c("-dbtype", "nucl", "-in", shQuote(ref_fasta), "-out", shQuote(db_prefix)),
  stdout = makeblastdb_stdout,
  stderr = makeblastdb_stderr
)
if (!identical(status, 0L)) fail(paste0("makeblastdb failed with exit code ", status, "."))

blast_path <- file.path(upstream_dir, "ref_blast.txt")
blast_stderr <- file.path(upstream_dir, "blastn.stderr.log")
status <- system2(
  blastn,
  args = c(
    "-db", shQuote(db_prefix),
    "-query", shQuote(seqs),
    "-evalue", "1e-20",
    "-max_target_seqs", "1000000",
    "-outfmt", "6",
    "-out", shQuote(blast_path),
    "-num_threads", as.character(threads)
  ),
  stdout = file.path(upstream_dir, "blastn.stdout.log"),
  stderr = blast_stderr
)
if (!identical(status, 0L)) fail(paste0("blastn failed with exit code ", status, "."))
if (!file.exists(blast_path) || file.info(blast_path)$size == 0) {
  fail("Tax4Fun2 BLAST found no reference match; no prediction is possible.")
}

blast <- read.delim(blast_path, header = FALSE, stringsAsFactors = FALSE)
if (ncol(blast) < 3) fail("Tax4Fun2 BLAST output is malformed.")
blast <- blast[!duplicated(blast[[1]]), , drop = FALSE]
write.table(
  blast,
  blast_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
writeLines(c("RefBlast", database_mode, seqs, as.character(Sys.time())), file.path(upstream_dir, "logfile1.txt"))

Tax4Fun2::makeFunctionalPrediction(
  path_to_otu_table = otu_table_path,
  path_to_reference_data = reference_dir,
  path_to_temp_folder = upstream_dir,
  database_mode = database_mode,
  normalize_by_copy_number = TRUE,
  min_identity_to_reference = min_identity,
  normalize_pathways = normalize_pathways
)

functional_source <- file.path(upstream_dir, "functional_prediction.txt")
pathway_source <- file.path(upstream_dir, "pathway_prediction.txt")
required_file(functional_source, "Tax4Fun2 functional prediction")
required_file(pathway_source, "Tax4Fun2 pathway prediction")
functional <- read.delim(functional_source, check.names = FALSE, stringsAsFactors = FALSE)
pathways <- read.delim(pathway_source, check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(functional) == 0 || names(functional)[1] != "KO") fail("Invalid Tax4Fun2 functional output.")
if (nrow(pathways) == 0 || names(pathways)[1] != "pathway") fail("Invalid Tax4Fun2 pathway output.")
write.table(
  functional,
  file.path(output_dir, "functional_prediction.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  pathways,
  file.path(output_dir, "pathway_prediction.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

identity_percent <- min_identity * 100
passing <- blast[blast[[3]] >= identity_percent, , drop = FALSE]
matched_ids <- unique(as.character(passing[[1]]))
matched <- feature_ids %in% matched_ids
positive <- abundance > 0
coverage <- data.frame(
  sample = sample_names,
  feature_fraction_used = colSums(positive[matched, , drop = FALSE]) / colSums(positive),
  sequence_fraction_used = colSums(abundance[matched, , drop = FALSE]) / colSums(abundance),
  stringsAsFactors = FALSE
)
write.table(
  coverage,
  file.path(output_dir, "coverage.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

used_subjects <- unique(as.character(passing[[2]]))
used_profiles <- file.path(reference_dir, database_mode, paste0(used_subjects, ".tbl.gz"))
used_profiles <- used_profiles[file.exists(used_profiles)]
fingerprint_files <- unique(c(required_files, used_profiles))
normalized_reference_dir <- normalizePath(reference_dir)
fingerprints <- lapply(fingerprint_files, function(path) {
  info <- file.info(path)
  normalized_path <- normalizePath(path)
  list(
    relative_path = substring(normalized_path, nchar(normalized_reference_dir) + 2),
    size = unname(info$size),
    md5 = unname(tools::md5sum(path))
  )
})

manifest <- list(
  schema_version = "microsuite-tax4fun2.v1",
  tax4fun2_version = as.character(utils::packageVersion("Tax4Fun2")),
  database_mode = database_mode,
  database_fingerprint = fingerprints,
  parameters = list(
    min_identity = min_identity,
    normalize_by_copy_number = TRUE,
    normalize_pathways = normalize_pathways,
    threads = threads
  ),
  input = list(features = length(feature_ids), samples = length(sample_names)),
  matched_features = length(intersect(feature_ids, matched_ids)),
  outputs = list(
    functions = "functional_prediction.tsv",
    pathways = "pathway_prediction.tsv",
    coverage = "coverage.tsv"
  )
)
jsonlite::write_json(
  manifest,
  file.path(output_dir, "tax4fun2_manifest.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
