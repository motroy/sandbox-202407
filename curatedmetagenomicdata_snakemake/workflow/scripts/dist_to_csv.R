d <- readRDS(snakemake@input[[1]])
write.csv(as.matrix(d), snakemake@output[[1]])
