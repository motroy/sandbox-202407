#!/usr/bin/env bash
# Link the 34-sample test subset (input/subset.tsv) out of genomes/ and stage ATB's AMRFinderPlus tables as the script's amr_dir.
# Run from ctxm_plasmid_test/ after 01_download_genomes.sh.
set -euo pipefail
mkdir -p subset_genomes amrfinder
while read -r s; do
  ln -sf "../genomes/${s}.fa.gz" "subset_genomes/${s}.fa.gz"
  cp "input/amrfinder_atb/${s}.amrfinder.tsv" "amrfinder/"
done < input/subset_biosamples.txt
ls subset_genomes | wc -l
