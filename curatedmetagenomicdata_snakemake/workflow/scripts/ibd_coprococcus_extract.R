# Notebook cell 25 (R part): Coprococcus genus- and species-level abundances
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages(library(SummarizedExperiment))
ibd_tse <- readRDS(snakemake@input[[1]])

genus_tse <- agglomerateByRank(ibd_tse, rank = "Genus")
gi <- grep("g__Coprococcus|Coprococcus", rownames(genus_tse), ignore.case = TRUE)
si <- grep("s__Coprococcus|Coprococcus", rownames(ibd_tse), ignore.case = TRUE)
if (length(gi) == 0 || length(si) == 0) stop("Could not find Coprococcus taxa in the dataset.")

out <- function(m, path) {
    df <- as.data.frame(t(m)); df$sample_id <- rownames(df); df$condition <- ibd_tse$condition
    write.csv(df, path, row.names = FALSE)
}
out(assay(genus_tse, "relative_abundance")[gi, , drop = FALSE], snakemake@output$genus)
out(assay(ibd_tse,  "relative_abundance")[si, , drop = FALSE], snakemake@output$species)
