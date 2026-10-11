#!/usr/bin/env bash
# Example: aggregate two Amira runs (E. coli ST131 SRR13116389; E. coli O103:H2 FWSEC0007 SRR7896269)
# and compare with AMRFinderPlus. Run from anywhere; needs only python3. Outputs go to ./output/.
set -euo pipefail
cd "$(dirname "$0")"
AGG=../../scripts/amira_aggregate.py
mkdir -p output

# 1) default: gene-family level, allele-count matrix, comparison with AMRFinderPlus
python3 -I $AGG ST131=input/ST131 FWSEC0007=input/FWSEC0007 \
  -o output/summary.tsv --long output/all_alleles.tsv --matrix output/matrix_alleles.tsv \
  --amrfinder ST131=input/ST131.amrfinderplus.tsv FWSEC0007=input/FWSEC0007.amrfinderplus.tsv \
  --compare output/compare_amrfinderplus.tsv

# 2) estimated copies per gene instead of raw allele counts
python3 -I $AGG ST131=input/ST131 FWSEC0007=input/FWSEC0007 \
  -o /dev/null --matrix output/matrix_copies.tsv --matrix-value copies

# 3) drop weak alleles (relative read depth < 0.7)
python3 -I $AGG ST131=input/ST131 FWSEC0007=input/FWSEC0007 \
  -o output/summary_min_depth_0.7.tsv --min-depth 0.7
