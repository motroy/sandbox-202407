# Download one curatedMetagenomicData study and write standardised tables.
source(file.path(snakemake@scriptdir, "common.R"))
suppressPackageStartupMessages(library(curatedMetagenomicData))
dir.create(dirname(snakemake@output$taxa), recursive = TRUE, showWarnings = FALSE)
ExperimentHub::setExperimentHubOption("CACHE", normalizePath(snakemake@config$experimenthub_cache))
study <- snakemake@wildcards$study

# Pick the most recent snapshot that provides the data type, then download only that one.
fetch <- function(type) {
    avail <- tryCatch(curatedMetagenomicData(paste0(study, ".", type), dryrun = TRUE), error = function(e) character())
    if (length(avail) == 0) return(NULL)
    latest <- sort(avail, decreasing = TRUE)[1]
    message(study, " ", type, ": using ", latest)
    suppressMessages(curatedMetagenomicData(latest, dryrun = FALSE))[[1]]
}

tax <- fetch("relative_abundance")
if (is.null(tax)) stop("Study '", study, "' not found in curatedMetagenomicData (relative_abundance)")
write_matrix(assay(tax, 1), snakemake@output$taxa)

meta <- as.data.frame(colData(tax))
meta <- data.frame(sample_id = rownames(meta), dataset = study, meta, check.names = FALSE, stringsAsFactors = FALSE)
write.table(meta, snakemake@output$meta, sep = "\t", quote = FALSE, row.names = FALSE)

path <- fetch("pathway_abundance")
if (is.null(path)) {
    message(study, ": no pathway_abundance available; writing an empty pathways table")
    write.table(data.frame(feature = character()), snakemake@output$pathways, sep = "\t", quote = FALSE, row.names = FALSE)
} else {
    m <- assay(path, 1)
    m <- m[!grepl("|", rownames(m), fixed = TRUE), , drop = FALSE]   # unstratified pathways only
    write_matrix(m, snakemake@output$pathways)
}
