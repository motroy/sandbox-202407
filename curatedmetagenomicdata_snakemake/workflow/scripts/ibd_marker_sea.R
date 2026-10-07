# Notebook cell 22-23: marker-clade enrichment (fgsea) on marker_abundance profiles
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages({library(SummarizedExperiment); library(fgsea)})
set_cache(snakemake@config$experimenthub_cache)
cfg <- snakemake@config$ibd
studies <- unlist(cfg$studies)

lst <- lapply(studies, load_cmd, type = "marker_abundance")
lst <- lst[!vapply(lst, is.null, logical(1))]
ibd_markers_tse <- mergeSEs(lst, assay.type = "relative_abundance", missing_values = 0, join = "inner")

colData(ibd_markers_tse)$condition <- ifelse(grepl("control|healthy", ibd_markers_tse$study_condition, ignore.case = TRUE), "Healthy", "IBD")
ibd_markers_tse <- ibd_markers_tse[, ibd_markers_tse$body_site == cfg$body_site & !is.na(ibd_markers_tse$condition)]
message("Marker TSE: ", nrow(ibd_markers_tse), " markers x ", ncol(ibd_markers_tse), " samples")

ibd_markers_tse <- transformAssay(ibd_markers_tse, assay.type = "relative_abundance", method = "log10",
                                  pseudocount = cfg$log_pseudocount_sea, name = "log_abundance")
stats <- apply(assay(ibd_markers_tse, "log_abundance"), 1, function(x) {
    fit <- tryCatch(lm(x ~ ibd_markers_tse$condition), error = function(e) NULL)
    if (is.null(fit)) return(0)
    summary(fit)$coefficients[2, "t value"]
})
stats <- stats[is.finite(stats)]

tax_table <- as.data.frame(rowData(ibd_markers_tse))
group_col <- if ("species" %in% colnames(tax_table)) tax_table$species else rownames(tax_table)
clade_sets <- split(rownames(tax_table), group_col)
message("Marker rowData columns: ", paste(colnames(tax_table), collapse = ", "), "; clade sets: ", length(clade_sets))

set.seed(snakemake@config$seed)
res <- fgsea(pathways = clade_sets, stats = stats,
             minSize = cfg$marker_sea$min_size, maxSize = cfg$marker_sea$max_size)
empty <- data.frame(pathway = character(), pval = numeric(), padj = numeric(), ES = numeric(), NES = numeric(), size = numeric())
if (!is.null(res) && nrow(res) > 0) {
    res <- as.data.frame(res[order(res$padj), ]); res$leadingEdge <- NULL
} else res <- empty
write.csv(res, snakemake@output$all, row.names = FALSE)
sig <- res[!is.na(res$padj) & res$padj < cfg$fdr_threshold, ]
write.csv(sig, snakemake@output$sig, row.names = FALSE)
print(head(res, 15))
message("Significant marker clades (FDR < ", cfg$fdr_threshold, "): ", nrow(sig))
