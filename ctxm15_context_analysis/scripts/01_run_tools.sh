#!/usr/bin/env bash
# Per-assembly tools consumed by 02/05. Usage: 01_run_tools.sh <sample> <assembly.fa>
# Needs: amrfinder (AMRFinderPlus, `amrfinder -u` done), mob_recon, mlst, ClermonTyping, ectyper, FimTyper.
set -euo pipefail
s=$1; fa=$2; t=${THREADS:-4}
mkdir -p amrfinder mobrecon typing
amrfinder -n "$fa" --plus --organism Escherichia --threads "$t" -o "amrfinder/$s.amrfinder.tsv"
mob_recon -i "$fa" -o "mobrecon/$s" -n "$t" --force        # contig_report.txt (chromosome/plasmid) + mge.report.txt
mlst "$fa" > "typing/$s.mlst.tsv"
ectyper -i "$fa" -o "typing/$s.ectyper" --cores "$t" && cp "typing/$s.ectyper/output.tsv" "typing/$s.ectyper.tsv"
# ClermonTyping and FimTyper: install-specific entry points; adapt to your setup. Expected outputs:
#   typing/$s.clermont.txt  (phylogroup table)   typing/$s.fimtyper.tsv  (FimTyper results_tab.tsv)
# clermonTyping.sh --fasta "$fa" --name "$s" --threshold 1   (then copy <s>_phylogroups.txt)
# fimtyper.py -i "$fa" -o typing/$s.fimtyper -d <fimtyper_db> -k 95.00 -l 0.60
