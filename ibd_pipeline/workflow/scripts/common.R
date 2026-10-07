# Shared helpers for R steps
suppressPackageStartupMessages(library(SummarizedExperiment))
if (exists("snakemake") && length(snakemake@log) > 0) {
    .log_con <- file(snakemake@log[[1]], open = "wt")
    sink(.log_con, type = "output", split = TRUE)
    sink(.log_con, type = "message")
}
read_matrix <- function(path) {
    m <- read.delim(path, check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)
    as.matrix(m)
}
write_matrix <- function(m, path) {
    out <- data.frame(feature = rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE)
    write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}
