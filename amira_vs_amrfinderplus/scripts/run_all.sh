#!/usr/bin/env bash
# Reproduce the Amira (SRR13116389 reads) vs AMRFinderPlus (GCF_025392675.1 assembly) comparison.
# Run from the amira_vs_amrfinderplus/ directory. Needs pixi (https://pixi.sh); ~1.5 GB download, ~25 min on 4 cores.
set -euo pipefail
mkdir -p bin data results/amira results/amrfinderplus results/summary

pixi install --all

# Pandora 0.12.0 (Amira needs this version; not on conda)
curl -fsSL -o bin/pandora https://github.com/iqbal-lab-org/pandora/releases/download/0.12.0-alpha.0/pandora-linux-precompiled-v0.12.0-alpha.0
chmod +x bin/pandora

# E. coli panRG (figshare file 54318899; use the ndownloader.figshare.com host -- www.figshare.com is bot-challenged)
curl -fsSL -o data/Escherichia.coli.panidx.zip https://ndownloader.figshare.com/files/54318899

# Assembly (Illumina/SPAdes, E. coli ST131)
curl -fsSL -o data/GCF_025392675.1.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/025/392/675/GCF_025392675.1_ASM2539267v1/GCF_025392675.1_ASM2539267v1_genomic.fna.gz
gzip -dc data/GCF_025392675.1.fna.gz > data/GCF_025392675.1.fna

# Reads: first 100,000 of 800,000 ONT reads (~683 Mb, ~130x) -- the full run is 5.7 GB compressed
curl -fsSL https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR131/089/SRR13116389/SRR13116389_1.fastq.gz | gzip -dc | head -n 400000 | gzip > data/SRR13116389_100k.fastq.gz || true

# AMRFinderPlus 4.2.7 + latest database
pixi run -e amrfinder amrfinder -u
pixi run -e amrfinder amrfinder -n data/GCF_025392675.1.fna -O Escherichia --plus --threads 4 \
  -o results/amrfinderplus/GCF_025392675.1.amrfinderplus.tsv \
  --nucleotide_output results/amrfinderplus/GCF_025392675.1.amr_genes.fna

# Amira 0.11.0
pixi run -e amira amira --reads data/SRR13116389_100k.fastq.gz --output results/amira/SRR13116389 \
  --species Escherichia_coli --panRG-path data/Escherichia.coli.panidx.zip \
  --pandora-path "$PWD/bin/pandora" --cores 4

# Compare
pixi run python -I scripts/compare.py results/amira/SRR13116389/amira_results.tsv \
  results/amrfinderplus/GCF_025392675.1.amrfinderplus.tsv results/summary
