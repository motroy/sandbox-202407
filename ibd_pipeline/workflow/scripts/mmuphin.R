# Batch-effect correction with MMUPHin::adjust_batch (batch = configured column, covariate = condition)
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages(library(MMUPHin))
cfg <- snakemake@config
bc <- cfg$batch_correction
set.seed(cfg$seed)

mat  <- read_matrix(snakemake@input$matrix)
meta <- read.delim(snakemake@input$meta, stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$sample_id
meta <- meta[colnames(mat), ]

batch <- bc$batch_column
tab <- table(meta[[batch]])
message("Batch sizes (", batch, "): ", paste(names(tab), tab, sep = "=", collapse = ", "))
if (length(tab) < 2) stop("MMUPHin needs at least 2 batches in '", batch, "'; use batch_correction.mode: none")
small <- names(tab)[tab < bc$min_batch_size]
if (length(small) > 0) stop("Batches smaller than batch_correction.min_batch_size (", bc$min_batch_size, "): ",
                            paste(small, collapse = ", "), ". Remove these datasets or lower min_batch_size.")
ct <- table(meta[[batch]], meta$condition)
if (any(ct == 0)) message("NOTE: some batches contain only one condition:\n", paste(capture.output(print(ct)), collapse = "\n"))

meta[[batch]]   <- factor(meta[[batch]])
meta$condition  <- factor(meta$condition, levels = unlist(cfg$condition_labels))

feat <- mat / 100                                  # MMUPHin expects proportions
keep <- rowMeans(feat > 0) >= bc$min_prevalence
message("Features kept (prevalence >= ", bc$min_prevalence, "): ", sum(keep), " of ", length(keep))
feat <- feat[keep, , drop = FALSE]

fit <- adjust_batch(feature_abd = feat, batch = batch, covariates = "condition", data = meta,
                    control = list(verbose = FALSE, diagnostic_plot = NULL))
adj <- fit$feature_abd_adj * 100
write_matrix(adj, snakemake@output$matrix)
write.csv(data.frame(feature = rownames(mat), kept = keep), snakemake@output$info, row.names = FALSE)
