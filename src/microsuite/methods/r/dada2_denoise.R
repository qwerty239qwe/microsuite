suppressPackageStartupMessages({
  library(dada2)
})

args <- commandArgs(trailingOnly = TRUE)

value_after <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) {
    return(default)
  }
  args[[idx + 1]]
}

has_flag <- function(flag) {
  flag %in% args
}

stem_without_fastq_suffix <- function(path) {
  name <- basename(path)
  sub("\\.(fastq|fq)(\\.gz)?$", "", name, ignore.case = TRUE)
}

read_suffix_pattern <- function(read) {
  direction <- if (read == 1) "forward" else "reverse"
  paste0(
    "([._-]R", read, "([._-]001)?",
    "|[._-]read", read, "([._-]001)?",
    "|[._-]", read, "([._-]001)?",
    "|[._-]?", direction,
    ")$"
  )
}

is_read <- function(path, read) {
  grepl(read_suffix_pattern(read), stem_without_fastq_suffix(path), ignore.case = TRUE)
}

sample_name <- function(path, read) {
  sub(read_suffix_pattern(read), "", stem_without_fastq_suffix(path), ignore.case = TRUE)
}

input_dir <- value_after("--input-dir")
output_table <- value_after("--output-table")
output_rep_seqs <- value_after("--output-rep-seqs")
output_stats <- value_after("--output-stats")
output_plot_dir <- value_after("--output-plot-dir")
threads <- as.integer(value_after("--threads", "1"))
paired <- has_flag("--paired")
max_n <- as.integer(value_after("--max-n", "0"))
trunc_q <- as.integer(value_after("--trunc-q", "2"))
rm_phix <- if (has_flag("--rm-phix")) TRUE else if (has_flag("--no-rm-phix")) FALSE else TRUE
pooling_method <- value_after("--pooling-method", "independent")
pool <- if (pooling_method == "pseudo") "pseudo" else FALSE
chimera_method <- value_after("--chimera-method", "consensus")
min_fold_parent_over_abundance <- as.numeric(value_after("--min-fold-parent-over-abundance", "1.0"))
allow_one_off <- has_flag("--allow-one-off")
n_reads_learn <- as.integer(value_after("--n-reads-learn", "1000000"))

params_out <- value_after("--params-out")

resolved_trim_left_f <- as.integer(value_after("--trim-left-f", "0"))
resolved_trim_left_r <- as.integer(value_after("--trim-left-r", "0"))
resolved_trunc_len_f <- as.integer(value_after("--trunc-len-f", "0"))
resolved_trunc_len_r <- as.integer(value_after("--trunc-len-r", "0"))
resolved_max_ee_f <- as.numeric(value_after("--max-ee-f", "2"))
resolved_max_ee_r <- as.numeric(value_after("--max-ee-r", "2"))
resolved_min_overlap <- as.integer(value_after("--min-overlap", "12"))
resolved_max_merge_mismatch <- as.integer(value_after("--max-merge-mismatch", "0"))
resolved_trim_overhang <- has_flag("--trim-overhang")
resolved_trim_left <- as.integer(value_after("--trim-left", "0"))
resolved_trunc_len <- as.integer(value_after("--trunc-len", "0"))
resolved_max_ee <- as.numeric(value_after("--max-ee", "2"))

if (is.null(input_dir) || is.null(output_table) || is.null(output_rep_seqs) || is.null(output_stats)) {
  stop("Missing required --input-dir, --output-table, --output-rep-seqs, or --output-stats.")
}
if (!pooling_method %in% c("independent", "pseudo")) {
  stop("--pooling-method must be independent or pseudo.")
}
if (!chimera_method %in% c("consensus", "none")) {
  stop("--chimera-method must be consensus or none.")
}

fastqs <- sort(list.files(input_dir, pattern = "\\.(fastq|fq)(\\.gz)?$", full.names = TRUE))
if (length(fastqs) == 0) {
  stop("No FASTQ files found in input directory.")
}

dir.create(dirname(output_table), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_rep_seqs), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output_stats), recursive = TRUE, showWarnings = FALSE)
if (!is.null(output_plot_dir)) {
  dir.create(output_plot_dir, recursive = TRUE, showWarnings = FALSE)
}

getN <- function(x) sum(getUniques(x))

write_error_plot <- function(err, output) {
  png(output, width = 1400, height = 1000, res = 150)
  on.exit(dev.off(), add = TRUE)
  print(plotErrors(err, nominalQ = TRUE))
}

write_retention_plot <- function(track, output) {
  png(output, width = 1600, height = 1000, res = 150)
  on.exit(dev.off(), add = TRUE)
  values <- as.matrix(track)
  input <- pmax(values[, "input"], 1)
  proportions <- sweep(values, 1, input, "/")
  matplot(
    t(proportions),
    type = "l",
    lty = 1,
    col = grDevices::hcl.colors(nrow(proportions), "Dark 3"),
    xaxt = "n",
    ylim = c(0, 1),
    xlab = "DADA2 step",
    ylab = "Fraction of input reads",
    main = "DADA2 read retention"
  )
  axis(1, at = seq_len(ncol(proportions)), labels = colnames(proportions), las = 2)
  grid()
}

# dada() and mergePairs() return a bare per-sample object for a single sample and
# a named list for multiple samples. Wrap the single-sample case so the read
# tallies below work either way (otherwise vapply() iterates the object's
# internal slots and getUniques() errors on them).
count_reads <- function(x) {
  if (is.data.frame(x) || inherits(x, "dada") || inherits(x, "derep")) {
    return(getN(x))
  }
  vapply(x, getN, numeric(1))
}

json_scalar <- function(v) {
  if (is.null(v) || (length(v) == 1 && is.na(v))) return("null")
  if (is.logical(v)) return(if (isTRUE(v)) "true" else "false")
  if (is.numeric(v)) return(format(v, scientific = FALSE, trim = TRUE))
  paste0("\"", gsub("\"", "\\\\\"", as.character(v)), "\"")
}

write_params_json <- function(path, params) {
  parts <- vapply(
    names(params),
    function(k) paste0("  \"", k, "\": ", json_scalar(params[[k]])),
    character(1)
  )
  writeLines(c("{", paste(parts, collapse = ",\n"), "}"), path)
}

if (paired) {
  fnFs <- fastqs[vapply(fastqs, is_read, logical(1), read = 1)]
  fnRs <- fastqs[vapply(fastqs, is_read, logical(1), read = 2)]
  sampleFs <- sample_name(fnFs, 1)
  sampleRs <- sample_name(fnRs, 2)
  if (length(fnFs) == 0 || length(fnFs) != length(fnRs) || !identical(sampleFs, sampleRs)) {
    stop("Paired DADA2 mode requires matching R1/R2 FASTQ files with consistent sample names.")
  }
  filtFs <- file.path(tempdir(), paste0(sampleFs, ".R1.filtered.fastq.gz"))
  filtRs <- file.path(tempdir(), paste0(sampleFs, ".R2.filtered.fastq.gz"))
  out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    trimLeft = c(resolved_trim_left_f, resolved_trim_left_r),
    truncLen = c(resolved_trunc_len_f, resolved_trunc_len_r),
    maxEE = c(resolved_max_ee_f, resolved_max_ee_r),
    truncQ = trunc_q,
    maxN = max_n,
    rm.phix = rm_phix,
    multithread = threads
  )
  errF <- learnErrors(filtFs, nbases = n_reads_learn, multithread = threads)
  errR <- learnErrors(filtRs, nbases = n_reads_learn, multithread = threads)
  if (!is.null(output_plot_dir)) {
    write_error_plot(errF, file.path(output_plot_dir, "error_rates_forward.png"))
    write_error_plot(errR, file.path(output_plot_dir, "error_rates_reverse.png"))
  }
  dadaFs <- dada(filtFs, err = errF, pool = pool, multithread = threads)
  dadaRs <- dada(filtRs, err = errR, pool = pool, multithread = threads)
  names(dadaFs) <- sampleFs
  names(dadaRs) <- sampleFs
  mergers <- mergePairs(
    dadaFs, filtFs, dadaRs, filtRs,
    minOverlap = resolved_min_overlap,
    maxMismatch = resolved_max_merge_mismatch,
    trimOverhang = resolved_trim_overhang
  )
  names(mergers) <- sampleFs
  seqtab <- makeSequenceTable(mergers)
  seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method = chimera_method,
    minFoldParentOverAbundance = min_fold_parent_over_abundance,
    allowOneOff = allow_one_off,
    multithread = threads
  )
  track <- data.frame(
    input = out[, "reads.in"],
    filtered = out[, "reads.out"],
    denoised_f = count_reads(dadaFs),
    denoised_r = count_reads(dadaRs),
    merged = count_reads(mergers),
    nonchim = rowSums(seqtab.nochim),
    row.names = sampleFs
  )
} else {
  samples <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(fastqs)))
  filt <- file.path(tempdir(), paste0(samples, ".filtered.fastq.gz"))
  out <- filterAndTrim(
    fastqs, filt,
    trimLeft = resolved_trim_left,
    truncLen = resolved_trunc_len,
    maxEE = resolved_max_ee,
    truncQ = trunc_q,
    maxN = max_n,
    rm.phix = rm_phix,
    multithread = threads
  )
  err <- learnErrors(filt, nbases = n_reads_learn, multithread = threads)
  if (!is.null(output_plot_dir)) {
    write_error_plot(err, file.path(output_plot_dir, "error_rates.png"))
  }
  dada_out <- dada(filt, err = err, pool = pool, multithread = threads)
  names(dada_out) <- samples
  seqtab <- makeSequenceTable(dada_out)
  seqtab.nochim <- removeBimeraDenovo(
    seqtab,
    method = chimera_method,
    minFoldParentOverAbundance = min_fold_parent_over_abundance,
    allowOneOff = allow_one_off,
    multithread = threads
  )
  track <- data.frame(
    input = out[, "reads.in"],
    filtered = out[, "reads.out"],
    denoised = count_reads(dada_out),
    nonchim = rowSums(seqtab.nochim),
    row.names = samples
  )
}

seqs <- colnames(seqtab.nochim)
ids <- paste0("ASV", seq_along(seqs))
asv_table <- t(seqtab.nochim)
rownames(asv_table) <- ids
write.table(asv_table, output_table, sep = "\t", quote = FALSE, col.names = NA)
writeLines(as.vector(rbind(paste0(">", ids), seqs)), output_rep_seqs)
write.table(track, output_stats, sep = "\t", quote = FALSE, col.names = NA)
if (!is.null(output_plot_dir)) {
  write_retention_plot(track, file.path(output_plot_dir, "read_retention.png"))
}

if (!is.null(params_out)) {
  write_params_json(params_out, list(
    mode = if (paired) "paired" else "single",
    trim_left_f = if (paired) resolved_trim_left_f else NA,
    trim_left_r = if (paired) resolved_trim_left_r else NA,
    trunc_len_f = if (paired) resolved_trunc_len_f else NA,
    trunc_len_r = if (paired) resolved_trunc_len_r else NA,
    max_ee_f = if (paired) resolved_max_ee_f else NA,
    max_ee_r = if (paired) resolved_max_ee_r else NA,
    min_overlap = if (paired) resolved_min_overlap else NA,
    max_merge_mismatch = if (paired) resolved_max_merge_mismatch else NA,
    trim_overhang = if (paired) resolved_trim_overhang else NA,
    trim_left = if (!paired) resolved_trim_left else NA,
    trunc_len = if (!paired) resolved_trunc_len else NA,
    max_ee = if (!paired) resolved_max_ee else NA,
    trunc_q = trunc_q,
    max_n = max_n,
    rm_phix = rm_phix,
    pooling_method = pooling_method,
    chimera_method = chimera_method,
    min_fold_parent_over_abundance = min_fold_parent_over_abundance,
    allow_one_off = allow_one_off,
    n_reads_learn = n_reads_learn,
    dada2_version = as.character(packageVersion("dada2")),
    r_version = R.version.string
  ))
}
