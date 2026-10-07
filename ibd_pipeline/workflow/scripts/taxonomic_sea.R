# fgsea over taxonomic groups, ranked by the differential-abundance t statistic (positive = case-enriched)
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages(library(fgsea))
cfg <- snakemake@config
da <- read.csv(snakemake@input[[1]], stringsAsFactors = FALSE)
stats <- setNames(da$t_value, da$feature)
stats <- stats[is.finite(stats)]

prefix <- c(Phylum = "p__", Class = "c__", Order = "o__", Family = "f__", Genus = "g__")
tax <- do.call(rbind, lapply(names(stats), function(f) {
    parts <- strsplit(f, "|", fixed = TRUE)[[1]]
    row <- vapply(prefix, function(p) { h <- parts[startsWith(parts, p)]; if (length(h)) sub(p, "", h[1], fixed = TRUE) else NA_character_ }, "")
    as.data.frame(as.list(row), stringsAsFactors = FALSE)
}))
rownames(tax) <- names(stats)

pathways <- list()
for (lvl in names(prefix)) {
    for (taxon in na.omit(unique(tax[[lvl]]))) {
        sp <- rownames(tax)[which(tax[[lvl]] == taxon)]
        if (length(sp) >= 1) pathways[[paste0(lvl, ":", taxon)]] <- sp
    }
}
set.seed(cfg$seed)
res <- fgsea(pathways = pathways, stats = stats, minSize = cfg$analysis$sea$min_size, maxSize = cfg$analysis$sea$max_size)
res <- as.data.frame(res[order(res$padj), ]); res$leadingEdge <- NULL
write.csv(res, snakemake@output$all, row.names = FALSE)
sig <- res[!is.na(res$padj) & res$padj < cfg$analysis$da$fdr_threshold, ]
write.csv(sig, snakemake@output$sig, row.names = FALSE)
message("Significant taxonomic groups: ", nrow(sig), " of ", nrow(res))
