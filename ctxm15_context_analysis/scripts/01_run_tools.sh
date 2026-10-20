#!/usr/bin/env bash
# Run the per-assembly tools consumed by 02_context_analysis.py.
# Usage: 01_run_tools.sh <sample> <assembly.fa> <is_references.fa>
# Needs: amrfinder, mob_recon, blastn, makeblastdb. is_references.fa holds ISEcp1, IS26 (e.g. from ISfinder).
set -euo pipefail
s=$1; fa=$2; isdb=$3
mkdir -p amrfinder mobrecon is_blast
amrfinder -n "$fa" --plus -o "amrfinder/$s.amrfinder.tsv"
mob_recon -i "$fa" -o "mobrecon/$s" -n "${THREADS:-4}" --force
[ -e "$isdb.nsq" ] || makeblastdb -in "$isdb" -dbtype nucl
# query = assembly contigs, subject db = IS references
blastn -query "$fa" -db "$isdb" -evalue 1e-20 \
  -outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue' > "is_blast/$s.tsv"
