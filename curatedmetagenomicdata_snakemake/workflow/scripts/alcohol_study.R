# Notebook cell 9: alcohol consumption & stool microbiota (adults, stool)
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages({library(SummarizedExperiment); library(stringr); library(ggplot2); library(uwot)})
set_cache(snakemake@config$experimenthub_cache)
set.seed(snakemake@config$seed)
min_age <- snakemake@config$alcohol$min_age

alcoholStudy <-
    filter(sampleMetadata, age >= min_age) |>
    filter(!is.na(alcohol)) |>
    filter(body_site == "stool") |>
    select(where(~ !all(is.na(.x)))) |>
    returnSamples("relative_abundance", rownames = "short")

colData(alcoholStudy) <-
    colData(alcoholStudy) |>
    as.data.frame() |>
    mutate(alcohol = str_replace_all(alcohol, "no", "No")) |>
    mutate(alcohol = str_replace_all(alcohol, "yes", "Yes")) |>
    DataFrame()

altExps(alcoholStudy) <- splitByRanks(alcoholStudy)

alcoholStudy <- addAlpha(alcoholStudy, assay.type = "relative_abundance", index = "shannon", name = "shannon")
p1 <- plotColData(alcoholStudy, x = "alcohol", y = "shannon", colour_by = "alcohol", shape_by = "alcohol") +
    labs(x = "Alcohol", y = "Alpha Diversity (H')") + theme(legend.position = "none")

alcoholStudy <- runMDS(alcoholStudy, FUN = vegan::vegdist, method = "bray",
                       exprs_values = "relative_abundance", altexp = "genus", name = "BrayCurtis")
p2 <- plotReducedDim(alcoholStudy, "BrayCurtis", colour_by = "alcohol", shape_by = "alcohol") +
    labs(x = "PCo 1", y = "PCo 2", title = "Bray-Curtis PCoA")

alcoholStudy <- runUMAP(alcoholStudy, exprs_values = "relative_abundance", altexp = "genus", name = "UMAP")
p3 <- plotReducedDim(alcoholStudy, "UMAP", colour_by = "alcohol", shape_by = "alcohol") +
    labs(x = "UMAP 1", y = "UMAP 2", title = "UMAP")

ggsave(snakemake@output$p_alpha, p1, width = 8, height = 6, dpi = 100)
ggsave(snakemake@output$p_bray,  p2, width = 8, height = 6, dpi = 100)
ggsave(snakemake@output$p_umap,  p3, width = 8, height = 6, dpi = 100)

# Tables for the python/seaborn re-plot (cell 11)
meta <- as.data.frame(colData(alcoholStudy))
meta$sample_id <- rownames(meta)
write.csv(meta[, c("sample_id", "study_name", "alcohol", "shannon")], snakemake@output$meta, row.names = FALSE)
bray <- as.data.frame(reducedDim(alcoholStudy, "BrayCurtis")); bray$sample_id <- rownames(meta)
umap <- as.data.frame(reducedDim(alcoholStudy, "UMAP"));       umap$sample_id <- rownames(meta)
write.csv(bray, snakemake@output$bray, row.names = FALSE)
write.csv(umap, snakemake@output$umap, row.names = FALSE)
