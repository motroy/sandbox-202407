# Bray-Curtis distance, PCoA, PERMANOVA (batch + condition, marginal) and beta-dispersion
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages(library(vegan))
cfg <- snakemake@config
set.seed(cfg$seed)
nperm <- cfg$analysis$permutations

mat  <- read_matrix(snakemake@input$matrix)
meta <- read.delim(snakemake@input$meta, stringsAsFactors = FALSE, check.names = FALSE)
rownames(meta) <- meta$sample_id
meta <- meta[colnames(mat), ]
batch_col <- cfg$batch_correction$batch_column
meta$batch <- factor(meta[[batch_col]])
meta$condition <- factor(meta$condition, levels = unlist(cfg$condition_labels))

d <- vegdist(t(mat), method = "bray")
pc <- cmdscale(d, k = 2, eig = TRUE)
ve <- 100 * pc$eig[1:2] / sum(abs(pc$eig))
coords <- data.frame(sample_id = rownames(pc$points), PCo1 = pc$points[, 1], PCo2 = pc$points[, 2],
                     condition = meta$condition, batch = meta$batch, var_PCo1 = ve[1], var_PCo2 = ve[2])
write.csv(coords, snakemake@output$coords, row.names = FALSE)

par <- snakemake@threads
tidy <- function(r, model) {
    data.frame(model = model, term = rownames(r), Df = r$Df, SumOfSqs = r$SumOfSqs, R2 = r$R2, F = r$F, p = r$`Pr(>F)`)
}
res <- tidy(adonis2(d ~ condition, data = meta, permutations = nperm, parallel = par), "condition_only")
if (nlevels(meta$batch) > 1) {
    r2 <- adonis2(d ~ batch + condition, data = meta, by = "margin", permutations = nperm, parallel = par)
    res <- rbind(tidy(r2, "batch+condition (marginal)"), res)
}
res <- res[!res$term %in% c("Residual", "Total"), ]
res$term[res$term == "batch"] <- batch_col
write.csv(res, snakemake@output$permanova, row.names = FALSE)
print(res)

bd <- betadisper(d, meta$condition)
pt <- permutest(bd, permutations = nperm)
disp <- data.frame(group = names(bd$group.distances), mean_distance_to_centroid = bd$group.distances,
                   F = pt$tab$F[1], p = pt$tab$`Pr(>F)`[1])
write.csv(disp, snakemake@output$dispersion, row.names = FALSE)
