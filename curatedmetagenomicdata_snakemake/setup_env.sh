#!/usr/bin/env bash
# One-off fix-ups after `pixi install` (idempotent). Run via: pixi run setup
#  1. curatedMetagenomicData 3.18.0 needs rbiom::unifrac, which was removed in rbiom >= 2,
#     so install rbiom 1.0.3 from the CRAN archive over the conda build.
#  2. The bioconda curatedMetagenomicData package downloads its payload in a post-link
#     script that pixi skips by default; run it by hand.
#  3. MMUPHin: the bioconda build only targets R 4.3, so install it from Bioconductor
#     (compiles several dependencies, ~20-30 min the first time).
set -euo pipefail
if ! Rscript -e 'quit(status = !requireNamespace("curatedMetagenomicData", quietly = TRUE))' 2>/dev/null; then
    Rscript -e 'if (packageVersion("rbiom") >= "2.0.0")
        install.packages("https://cran.r-project.org/src/contrib/Archive/rbiom/rbiom_1.0.3.tar.gz", repos = NULL, type = "source")'
    PREFIX="$CONDA_PREFIX" bash "$CONDA_PREFIX/bin/.bioconductor-curatedmetagenomicdata-post-link.sh"
fi
Rscript -e 'library(curatedMetagenomicData); cat("curatedMetagenomicData", as.character(packageVersion("curatedMetagenomicData")), "OK\n")'
if ! Rscript -e 'quit(status = !requireNamespace("MMUPHin", quietly = TRUE))' 2>/dev/null; then
    Rscript -e 'options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 1); BiocManager::install("MMUPHin", update = FALSE, ask = FALSE)'
fi
Rscript -e 'library(MMUPHin); cat("MMUPHin", as.character(packageVersion("MMUPHin")), "OK\n")'
