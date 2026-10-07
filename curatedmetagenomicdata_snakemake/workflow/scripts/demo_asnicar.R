# Notebook cells 4-5: AsnicarF_2017 demo -> PCA, Shannon, PERMANOVA, simple lm
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages(library(SummarizedExperiment))
set_cache(snakemake@config$experimenthub_cache)
set.seed(snakemake@config$seed)

study <- snakemake@config$demo$dataset
tse <- curatedMetagenomicData(paste0(study, ".relative_abundance"), dryrun = FALSE)[[1]]

# runPCA() needs a 'logcounts' assay, which curatedMetagenomicData does not ship.
tse <- logNormCounts(tse, assay.type = "relative_abundance")
tse <- runPCA(tse)

pca <- as.data.frame(reducedDim(tse, "PCA"))
pca$sample_id <- colnames(tse)
pca$subject_id <- colData(tse)$subject_id
pca$body_site <- colData(tse)$body_site
write.csv(pca, snakemake@output$pca, row.names = FALSE)

# Alpha diversity
tse <- addAlpha(tse, assay.type = "relative_abundance", index = "shannon", name = "shannon")
alpha <- as.data.frame(colData(tse)[, c("subject_id", "body_site", "shannon")])
alpha$sample_id <- rownames(alpha)
write.csv(alpha, snakemake@output$alpha, row.names = FALSE)

# Beta diversity: Bray-Curtis PERMANOVA by body site
dist_matrix <- vegdist(t(assay(tse, "relative_abundance")), method = "bray")
permanova <- adonis2(dist_matrix ~ body_site, data = as.data.frame(colData(tse)),
                     permutations = 999)

# Simple differential abundance example: first feature vs body site
da <- summary(lm(assay(tse, "logcounts")[1, ] ~ colData(tse)$body_site))

sink(snakemake@output$stats)
cat("--- Beta Diversity (PERMANOVA, Bray-Curtis ~ body_site) ---\n"); print(permanova)
cat("\n--- Differential Abundance Example (first feature, lm on logcounts) ---\n")
cat("Feature:", rownames(tse)[1], "\n"); print(da)
sink()
