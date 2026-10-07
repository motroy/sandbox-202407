# Notebook cell 20: taxonomic set enrichment (fgsea) of IBD vs Healthy
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages({library(SummarizedExperiment); library(fgsea)})
cfg <- snakemake@config$ibd
ibd_tse <- readRDS(snakemake@input[[1]])

ibd_tse <- transformAssay(ibd_tse, assay.type = "relative_abundance", method = "log10",
                          pseudocount = cfg$log_pseudocount_sea, name = "log_abundance")
rank_stats <- apply(assay(ibd_tse, "log_abundance"), 1, function(x) {
    fit <- tryCatch(lm(x ~ ibd_tse$condition), error = function(e) NULL)
    if (is.null(fit)) return(0)
    summary(fit)$coefficients[2, "t value"]   # positive = higher in IBD
})
rank_stats <- rank_stats[is.finite(rank_stats)]

tax_table <- as.data.frame(rowData(ibd_tse))
pathways <- list()
for (lvl in c("Phylum", "Class", "Order", "Family", "Genus")) {
    taxa <- unique(tax_table[[lvl]])
    taxa <- taxa[!is.na(taxa) & taxa != ""]
    for (taxon in taxa) {
        sp <- rownames(tax_table)[tax_table[[lvl]] == taxon]
        sp <- sp[sp %in% names(rank_stats)]
        if (length(sp) >= 1) pathways[[paste0(lvl, ":", taxon)]] <- sp
    }
}

set.seed(snakemake@config$seed)
res <- fgsea(pathways = pathways, stats = rank_stats,
             minSize = cfg$taxonomic_sea$min_size, maxSize = cfg$taxonomic_sea$max_size)
res <- as.data.frame(res[order(res$padj), ])
res$leadingEdge <- NULL
write.csv(res, snakemake@output$all, row.names = FALSE)
sig <- res[!is.na(res$padj) & res$padj < cfg$fdr_threshold, ]
write.csv(sig, snakemake@output$sig, row.names = FALSE)
print(head(res, 15))
message("Significant taxonomic groups (FDR < ", cfg$fdr_threshold, "): ", nrow(sig))
