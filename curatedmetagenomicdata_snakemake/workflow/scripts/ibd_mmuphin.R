# MMUPHin batch-effect correction (adjust_batch) of the IBD taxonomic profiles.
# batch = study_name, covariate of interest = condition (Healthy/IBD) so the disease
# signal is preserved while between-study effects are removed.
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages({library(SummarizedExperiment); library(MMUPHin)})
cfg <- snakemake@config$ibd
set.seed(snakemake@config$seed)
ibd_tse <- readRDS(snakemake@input[[1]])

# relative abundances are in percent; MMUPHin expects proportions
feat <- assay(ibd_tse, "relative_abundance") / 100
meta <- as.data.frame(colData(ibd_tse))
meta$study_name <- factor(meta$study_name)
meta$condition  <- factor(meta$condition, levels = c("Healthy", "IBD"))

# Features absent from (nearly) all samples break the per-batch fits
keep <- rowMeans(feat > 0) >= snakemake@config$mmuphin$min_prevalence
feat <- feat[keep, ]
message("Features kept (prevalence >= ", snakemake@config$mmuphin$min_prevalence, "): ", sum(keep))

fit <- adjust_batch(feature_abd = feat, batch = "study_name", covariates = "condition",
                    data = meta, control = list(verbose = FALSE, diagnostic_plot = NULL))
adj <- fit$feature_abd_adj

# distances + PCoA before/after
pc <- function(m) {
    d <- vegdist(t(m), method = "bray")
    list(d = d, pcoa = cmdscale(d, k = 2, eig = TRUE))
}
b <- pc(feat); a <- pc(adj)
ve <- function(p) round(100 * p$pcoa$eig[1:2] / sum(abs(p$pcoa$eig)), 1)
coords <- rbind(
    data.frame(stage = "Before", sample_id = colnames(feat), PCo1 = b$pcoa$points[, 1], PCo2 = b$pcoa$points[, 2]),
    data.frame(stage = "After",  sample_id = colnames(adj),  PCo1 = a$pcoa$points[, 1], PCo2 = a$pcoa$points[, 2]))
coords$study_name <- meta$study_name[match(coords$sample_id, rownames(meta))]
coords$condition  <- meta$condition[match(coords$sample_id, rownames(meta))]
coords$varexp1 <- ifelse(coords$stage == "Before", ve(b)[1], ve(a)[1])
coords$varexp2 <- ifelse(coords$stage == "Before", ve(b)[2], ve(a)[2])
write.csv(coords, snakemake@output$coords, row.names = FALSE)

# PERMANOVA: variance explained by study vs condition, before and after
perm <- function(d, stage) {
    r <- adonis2(d ~ study_name + condition, data = meta, by = "margin",
                 permutations = snakemake@config$mmuphin$permutations)
    data.frame(stage = stage, term = rownames(r), Df = r$Df, SumOfSqs = r$SumOfSqs, R2 = r$R2, F = r$F, p = r$`Pr(>F)`)
}
res <- rbind(perm(b$d, "Before"), perm(a$d, "After"))
res <- res[res$term %in% c("study_name", "condition"), ]
write.csv(res, snakemake@output$permanova, row.names = FALSE)
print(res)

write.csv(adj * 100, snakemake@output$adjusted)   # back to percent
