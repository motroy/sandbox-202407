#!/usr/bin/env bash
# One-off fix-ups after `pixi install` (idempotent). Run via: pixi run setup
#  1. curatedMetagenomicData 3.18.0 needs rbiom::unifrac, which was removed in rbiom >= 2,
#     so install rbiom 1.0.3 from the CRAN archive over the conda build.
#  2. The bioconda curatedMetagenomicData package downloads its payload in a post-link
#     script that pixi skips by default; run it by hand.
set -euo pipefail
if ! Rscript -e 'quit(status = !requireNamespace("curatedMetagenomicData", quietly = TRUE))' 2>/dev/null; then
    Rscript -e 'if (packageVersion("rbiom") >= "2.0.0")
        install.packages("https://cran.r-project.org/src/contrib/Archive/rbiom/rbiom_1.0.3.tar.gz", repos = NULL, type = "source")'
    PREFIX="$CONDA_PREFIX" bash "$CONDA_PREFIX/bin/.bioconductor-curatedmetagenomicdata-post-link.sh"
fi
Rscript -e 'library(curatedMetagenomicData); cat("curatedMetagenomicData", as.character(packageVersion("curatedMetagenomicData")), "OK\n")'
