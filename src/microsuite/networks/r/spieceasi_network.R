args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("usage: spieceasi_network.R <sample-feature-table.tsv> <output.tsv> <method> <lambda-min-ratio> <nlambda>", call. = FALSE)
}

table_path <- args[[1]]
output <- args[[2]]
method <- args[[3]]
lambda_min_ratio <- as.numeric(args[[4]])
nlambda <- as.integer(args[[5]])

if (!requireNamespace("SpiecEasi", quietly = TRUE)) {
  stop("SpiecEasi R package is required.", call. = FALSE)
}

counts <- read.table(
  table_path,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  comment.char = "",
  quote = ""
)

fit <- SpiecEasi::spiec.easi(
  as.matrix(counts),
  method = method,
  lambda.min.ratio = lambda_min_ratio,
  nlambda = nlambda
)
adj <- as.matrix(SpiecEasi::getRefit(fit))
features <- colnames(counts)

edges <- data.frame(
  source = character(),
  target = character(),
  weight = numeric(),
  abs_weight = numeric(),
  stringsAsFactors = FALSE
)
for (i in seq_len(ncol(adj) - 1)) {
  for (j in seq.int(i + 1, ncol(adj))) {
    weight <- adj[i, j]
    if (!is.na(weight) && weight != 0) {
      edges <- rbind(edges, data.frame(
        source = features[[i]],
        target = features[[j]],
        weight = weight,
        abs_weight = abs(weight),
        stringsAsFactors = FALSE
      ))
    }
  }
}

write.table(edges, output, sep = "\t", row.names = FALSE, quote = FALSE)
