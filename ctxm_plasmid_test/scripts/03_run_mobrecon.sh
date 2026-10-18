#!/usr/bin/env bash
# Run MOB-recon on the subset assemblies, all CTX-M-15 carriers first (4 in parallel, 1 thread each) -> mobrecon/<SAMN>/
# Run from ctxm_plasmid_test/ with the `ctxm` env (environment.yml + mob_suite) on PATH.
set -uo pipefail
mkdir -p mobrecon fasta_tmp
run() {
  s=$1
  [ -s mobrecon/$s/contig_report.txt ] && exit 0
  zcat subset_genomes/$s.fa.gz > fasta_tmp/$s.fasta
  mob_recon -i fasta_tmp/$s.fasta -o mobrecon/$s -n 1 --force > mobrecon_$s.log 2>&1
  rm -f fasta_tmp/$s.fasta
}
export -f run
xargs -a input/mobrecon_order.txt -P 4 -I{} bash -c 'run {}'
