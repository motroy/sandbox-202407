# Notebook cells 13-16: IBD meta-analysis data assembly
#   * availability of functional data types
#   * load + merge 5 studies (taxonomic profiles), define Healthy/IBD
#   * Shannon + Faith PD, Bray-Curtis distance + PCoA
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages(library(SummarizedExperiment))
set_cache(snakemake@config$experimenthub_cache)
set.seed(snakemake@config$seed)
cfg <- snakemake@config$ibd
studies <- unlist(cfg$studies)

avail <- data.frame(
    study_name = studies,
    pathways = check_type(studies, "pathabundance"),
    gene_families = check_type(studies, "genefamilies"),
    markers = check_type(studies, "marker_abundance"),
    row.names = NULL)
print(avail)
write.csv(avail, snakemake@output$availability, row.names = FALSE)

tse_list <- lapply(studies, load_cmd, type = "relative_abundance")
tse_list <- tse_list[!vapply(tse_list, is.null, logical(1))]
ibd_tse <- mergeSEs(tse_list, assay.type = "relative_abundance", missing_values = 0, join = "inner")

colData(ibd_tse)$condition <- ifelse(grepl("control|healthy", ibd_tse$study_condition, ignore.case = TRUE), "Healthy", "IBD")
ibd_tse <- ibd_tse[, ibd_tse$body_site == cfg$body_site & !is.na(ibd_tse$condition)]
ibd_tse <- ibd_tse[, colSums(assay(ibd_tse)) > 0]   # remove empty samples
message("Samples: ", ncol(ibd_tse), "  Features: ", nrow(ibd_tse))
print(table(ibd_tse$condition))

ibd_tse <- addAlpha(ibd_tse, assay.type = "relative_abundance", index = "shannon", name = "shannon_index")
ibd_tse <- addAlpha(ibd_tse, assay.type = "relative_abundance", index = "faith", name = "faiths_pd")

ibd_tse <- runMDS(ibd_tse, FUN = vegan::vegdist, method = "bray",
                  assay.type = "relative_abundance", name = "PCoA_Bray")
dist_matrix <- vegan::vegdist(t(assay(ibd_tse, "relative_abundance")), method = "bray")

saveRDS(ibd_tse, snakemake@output$tse)
write.csv(as.data.frame(colData(ibd_tse)), snakemake@output$meta)
saveRDS(dist_matrix, snakemake@output$dist)
write.csv(as.data.frame(reducedDim(ibd_tse, "PCoA_Bray")), snakemake@output$pcoa)
