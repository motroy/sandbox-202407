#!/usr/bin/env bash
# cgMLST distances for duplicate removal (chewBBACA v3.5.2 AlleleCall + cgmlst-dists v0.4.0).
# Usage: SCHEMA=<E. coli Achtman cgMLST schema dir, 2,513 core genes> 00b_cgmlst.sh <genomes_dir> <outdir>
# Download the schema from Chewie-NS with chewBBACA (see its docs); not fetched here.
set -euo pipefail
g=$1; o=$2
chewBBACA.py AlleleCall -i "$g" -g "$SCHEMA" -o "$o" --cpu "${THREADS:-8}"
cgmlst-dists "$o"/results_*/results_alleles.tsv > "$o/cgmlst_dists.tsv"
