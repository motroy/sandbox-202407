# Builds small *local-format* test data from real curatedMetagenomicData samples, so the
# local-input code paths (per-sample MetaPhlAn4 profiles, HUMAnN per-sample files, merged
# MetaPhlAn tables, metadata) can be tested without private data.
#   cohort_A (IjazUZ_2017): per-sample MetaPhlAn4-style profiles + HUMAnN pathabundance files
#   cohort_B (LiJ_2014)   : one merged MetaPhlAn table, no HUMAnN
# Run:  pixi run Rscript tests/make_fixtures.R
suppressPackageStartupMessages({library(curatedMetagenomicData); library(SummarizedExperiment)})
set.seed(1)
out <- "tests/fixtures"
get <- function(name) suppressMessages(curatedMetagenomicData(name, dryrun = FALSE))[[1]]
pick <- function(tse, n_ctl, n_case, ctl_label = "control", case_label = "IBD") {
    cd <- colData(tse); s <- cd$study_condition
    c(sample(colnames(tse)[s == ctl_label], n_ctl), sample(colnames(tse)[s == case_label], n_case))
}

# --- cohort A: per-sample MetaPhlAn4-format profiles
A <- get("2021-10-14.IjazUZ_2017.relative_abundance"); ids <- pick(A, 12, 12)
dir.create(file.path(out, "cohort_A/metaphlan"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "cohort_A/humann"), recursive = TRUE, showWarnings = FALSE)
profile_lines <- function(ab) {                       # ab: named numeric, names = full lineages
    ab <- ab[ab > 0]; parts <- strsplit(names(ab), "|", fixed = TRUE)
    lines <- character()
    for (r in 1:6) {                                  # higher ranks aggregated like MetaPhlAn does
        pre <- vapply(parts, function(p) paste(p[1:r], collapse = "|"), "")
        s <- tapply(ab, pre, sum); lines <- c(lines, sprintf("%s\t0\t%.5f\t", names(s), s))
    }
    sp <- sprintf("%s\t0\t%.5f\t", names(ab), ab)
    sgb <- sprintf("%s|t__SGB%d\t\t%.5f\t", names(ab), seq_along(ab), ab)   # SGB-level leaves (must be ignored)
    c(lines, sp, sgb)
}
for (id in ids) {
    writeLines(c("#mpa_vOct22_CHOCOPhlAnSGB_202212", "#test fixture derived from curatedMetagenomicData",
                 "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species", profile_lines(assay(A)[, id])),
               file.path(out, "cohort_A/metaphlan", paste0(id, "_profile.txt")))
}
P <- get("2021-03-31.IjazUZ_2017.pathway_abundance")
strat <- which(grepl("|", rownames(P), fixed = TRUE))[1:30]
rows <- c(which(!grepl("|", rownames(P), fixed = TRUE)), strat)
for (id in ids) {
    v <- assay(P)[rows, id]
    writeLines(c(paste0("# Pathway\t", id, "_Abundance"), sprintf("%s\t%.6f", rownames(P)[rows], v)),
               file.path(out, "cohort_A/humann", paste0(id, "_pathabundance.tsv")))
}
cdA <- as.data.frame(colData(A))[ids, ]
write.table(data.frame(sample_id = ids, diagnosis = ifelse(cdA$study_condition == "control", "healthy", "IBD"),
                       age = cdA$age, gender = cdA$gender),
            file.path(out, "cohort_A/metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# --- cohort B: one merged MetaPhlAn table (merge_metaphlan_tables.py layout)
B <- get("2021-03-31.LiJ_2014.relative_abundance"); idsB <- pick(B, 10, 14)
dir.create(file.path(out, "cohort_B"), recursive = TRUE, showWarnings = FALSE)
mB <- assay(B)[, idsB]; mB <- mB[rowSums(mB) > 0, ]
con <- file(file.path(out, "cohort_B/metaphlan_merged.tsv"), "w")
writeLines(c("#mpa_vOct22_CHOCOPhlAnSGB_202212", paste(c("clade_name", idsB), collapse = "\t")), con)
write.table(data.frame(rownames(mB), round(mB, 5)), con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)
cdB <- as.data.frame(colData(B))[idsB, ]
write.table(data.frame(SampleName = idsB, Group = ifelse(cdB$study_condition == "control", "Control", "Patient"),
                       age = cdB$age),
            file.path(out, "cohort_B/metadata.csv"), sep = ",", quote = FALSE, row.names = FALSE)
message("fixtures written to ", out)
