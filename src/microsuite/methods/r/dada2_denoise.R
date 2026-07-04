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

sample_name <- function(path, marker) {
  name <- basename(path)
  name <- sub("\\.(fastq|fq)(\\.gz)?$", "", name, ignore.case = TRUE)
  sub(marker, "", name, ignore.case = TRUE)
}

input_dir <- value_after("--input-dir")
output_table <- value_after("--output-table")
output_rep_seqs <- value_after("--output-rep-seqs")
output_stats <- value_after("--output-stats")
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

getN <- function(x) sum(getUniques(x))

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

if (paired) {
  fnFs <- fastqs[grepl("(_R1_|_1\\.|forward)", basename(fastqs), ignore.case = TRUE)]
  fnRs <- fastqs[grepl("(_R2_|_2\\.|reverse)", basename(fastqs), ignore.case = TRUE)]
  sampleFs <- sample_name(fnFs, "(_R1_|_1\\.|forward)")
  sampleRs <- sample_name(fnRs, "(_R2_|_2\\.|reverse)")
  if (length(fnFs) == 0 || length(fnFs) != length(fnRs) || !identical(sampleFs, sampleRs)) {
    stop("Paired DADA2 mode requires matching R1/R2 FASTQ files with consistent sample names.")
  }
  filtFs <- file.path(tempdir(), paste0(sampleFs, ".R1.filtered.fastq.gz"))
  filtRs <- file.path(tempdir(), paste0(sampleFs, ".R2.filtered.fastq.gz"))
  out <- filterAndTrim(
    fnFs, filtFs, fnRs, filtRs,
    trimLeft = c(as.integer(value_after("--trim-left-f", "0")), as.integer(value_after("--trim-left-r", "0"))),
    truncLen = c(as.integer(value_after("--trunc-len-f", "0")), as.integer(value_after("--trunc-len-r", "0"))),
    maxEE = c(as.numeric(value_after("--max-ee-f", "2")), as.numeric(value_after("--max-ee-r", "2"))),
    truncQ = trunc_q,
    maxN = max_n,
    rm.phix = rm_phix,
    multithread = threads
  )
  errF <- learnErrors(filtFs, nbases = n_reads_learn, multithread = threads)
  errR <- learnErrors(filtRs, nbases = n_reads_learn, multithread = threads)
  dadaFs <- dada(filtFs, err = errF, pool = pool, multithread = threads)
  dadaRs <- dada(filtRs, err = errR, pool = pool, multithread = threads)
  mergers <- mergePairs(
    dadaFs, filtFs, dadaRs, filtRs,
    minOverlap = as.integer(value_after("--min-overlap", "12")),
    maxMismatch = as.integer(value_after("--max-merge-mismatch", "0")),
    trimOverhang = has_flag("--trim-overhang")
  )
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
    trimLeft = as.integer(value_after("--trim-left", "0")),
    truncLen = as.integer(value_after("--trunc-len", "0")),
    maxEE = as.numeric(value_after("--max-ee", "2")),
    truncQ = trunc_q,
    maxN = max_n,
    rm.phix = rm_phix,
    multithread = threads
  )
  err <- learnErrors(filt, nbases = n_reads_learn, multithread = threads)
  dada_out <- dada(filt, err = err, pool = pool, multithread = threads)
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
