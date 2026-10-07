# Shared helpers for all R steps.
suppressPackageStartupMessages({
    library(curatedMetagenomicData)
    library(mia)
    library(scater)
    library(vegan)
    library(dplyr)
})

# Keep the ExperimentHub cache inside the project so runs are self-contained.
set_cache <- function(path) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    ExperimentHub::setExperimentHubOption("CACHE", normalizePath(path))
}

# Download one dataset from curatedMetagenomicData and give the first assay a
# uniform name (as done in the notebook). Returns NULL on failure.
load_cmd <- function(study, type, assay_name = "relative_abundance") {
    res <- tryCatch(
        curatedMetagenomicData(paste0(study, ".", type), dryrun = FALSE),
        error = function(e) { message("Failed to load ", study, ".", type, ": ", conditionMessage(e)); NULL })
    if (is.null(res)) return(NULL)
    tse <- res[[1]]
    assayNames(tse)[1] <- assay_name
    tse
}

# TRUE/FALSE availability of a data type for a study (dry run, no download)
check_type <- function(studies, type) {
    vapply(studies, function(s) {
        res <- tryCatch(curatedMetagenomicData(paste0(s, ".", type), dryrun = TRUE),
                        error = function(e) NULL)
        !is.null(res) && length(res) > 0
    }, logical(1))
}

# Send stdout + messages/warnings to the rule's log file (script: jobs don't do this by default)
if (exists("snakemake") && length(snakemake@log) > 0) {
    .log_con <- file(snakemake@log[[1]], open = "wt")
    sink(.log_con, type = "output", split = TRUE)
    sink(.log_con, type = "message")
}
