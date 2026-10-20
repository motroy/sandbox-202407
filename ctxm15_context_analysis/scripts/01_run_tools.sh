#!/usr/bin/env bash
# Per-assembly tools consumed by 02/05. Usage: 01_run_tools.sh <sample> <assembly.fa>
# Needs: amrfinder (AMRFinderPlus, `amrfinder -u` done), mob_recon, mlst, clermonTyping.sh, ectyper, FimTyper (see below).
set -euo pipefail
s=$1; fa=$2; t=${THREADS:-4}
mkdir -p amrfinder mobrecon typing
amrfinder -n "$fa" --plus --organism Escherichia --threads "$t" -o "amrfinder/$s.amrfinder.tsv"
mob_recon -i "$fa" -o "mobrecon/$s" -n "$t" --force        # contig_report.txt (chromosome/plasmid) + mge.report.txt
mlst "$fa" > "typing/$s.mlst.tsv"
ectyper -i "$fa" -o "typing/$s.ectyper" --cores "$t" && cp "typing/$s.ectyper/output.tsv" "typing/$s.ectyper.tsv"
# ClermonTyping (bioconda `clermontyping`; its R report step may print an error, the phylogroup table is still written)
fa_abs=$(realpath "$fa")
(cd typing && clermonTyping.sh --fasta "$fa_abs" --name "$s" --threshold 1 >/dev/null 2>&1; cp "$s/${s}_phylogroups.txt" "$s.clermont.txt")
# FimTyper (Perl; needs bioperl, Try::Tiny::Retry, blast). FIMTYPER_DIR = clone of bitbucket.org/genomicepidemiology/fimtyper,
# FIMTYPER_DB = clone of .../fimtyper_db. -b takes a BLAST prefix with bin/blastn (the conda env).
perl "${FIMTYPER_DIR:?set FIMTYPER_DIR}/fimtyper.pl" -d "${FIMTYPER_DB:?set FIMTYPER_DB}" -b "${CONDA_PREFIX:?activate the env}" \
  -i "$fa_abs" -o "typing/$s.fimtyper" -k 95 -l 0.6 >/dev/null 2>&1
cp "typing/$s.fimtyper/results_tab.txt" "typing/$s.fimtyper.tsv"
