#!/usr/bin/env bash
# Per-assembly tools consumed by 02/05. Usage: 01_run_tools.sh <sample> <assembly.fa>
# Needs: amrfinder (AMRFinderPlus, `amrfinder -u` done), mob_recon, mlst, clermonTyping.sh, ectyper, FimTyper (see below).
set -euo pipefail
s=$1; fa=$2; t=${THREADS:-4}
mkdir -p amrfinder mobrecon typing resfinder platon isfinder
# ResFinder 4 (acquired genes, 80% identity/coverage; RESFINDER_DB = clone of bitbucket.org/genomicepidemiology/resfinder_db) and Platon
# (PLATON_DB = db.tar.gz from zenodo.org/records/4066768) are the default ARG / chromosome-plasmid callers (see config.toml).
run_resfinder.py -ifa "$fa" -o "resfinder/$s" -db_res "${RESFINDER_DB:?set RESFINDER_DB}" -acq -l 0.8 -t 0.8 -b "$(command -v blastn)" >/dev/null
platon --db "${PLATON_DB:?set PLATON_DB}" --threads "$t" --output "platon/$s" "$fa" >/dev/null 2>&1 || true   # no table is written when no contig passes Platon's filters
# IS elements: ABRicate vs ISfinder (custom db, see README); ISEc9 = ISEcp1 is handled in step 02
abricate --db isfinder --minid 80 --mincov 80 --quiet "$fa" > "isfinder/$s.tsv"
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
