#!/usr/bin/env bash
# Phylogeny of the genome collection: Mashtree (fast, assembly-based), Panaroo pangenome -> core-gene alignment -> IQ-TREE.
# Usage: 08_phylogeny.sh <genomes_dir> <outdir>      (genomes_dir holds <id>.fa and <id>.gbff)
# Needs: mashtree, panaroo (+ mafft), iqtree3 (or iqtree). Settings (assumptions, not from the paper) via env vars below.
set -euo pipefail
g=$1; o=$2; t=${THREADS:-4}
CORE=${CORE_THRESHOLD:-0.98}        # Panaroo: gene present in >= this fraction of genomes counts as core
CLEAN=${CLEAN_MODE:-strict}         # Panaroo graph cleaning
BOOT=${BOOTSTRAPS:-1000}            # IQ-TREE ultrafast bootstrap replicates
mkdir -p "$o"
IQ=$(command -v iqtree3 || command -v iqtree)

# 1. Mashtree: tip labels = file prefixes
mashtree --numcpus "$t" --outtree "$o/mashtree.dnd" "$g"/*.fa

# 2. Panaroo on the NCBI GBFFs (converted internally, so no re-annotation). Core alignment with MAFFT.
# Panaroo converts GenBank files only when -i is ONE text file listing the paths (passing *.gbff directly makes it
# parse them as GFF3 and fail). Paths must not contain spaces. NCBI annotations contain frameshifted CDS, which
# Panaroo rejects unless --remove-invalid-genes drops them.
ls "$g"/*.gbff > "$o/panaroo_input.txt"
panaroo -i "$o/panaroo_input.txt" -o "$o/panaroo" --clean-mode "$CLEAN" --remove-invalid-genes -a core --aligner mafft --core_threshold "$CORE" -t "$t"

# 3. IQ-TREE: ModelFinder + ultrafast bootstrap on the core-gene alignment.
# Invariant-site-aware: Panaroo's alignment is not SNP-only, so no ASC correction.
"$IQ" -s "$o/panaroo/core_gene_alignment.aln" -m MFP -B "$BOOT" -bnni -T AUTO --threads-max "$t" --prefix "$o/iqtree" -redo
echo "trees: $o/mashtree.dnd  $o/iqtree.treefile  (pangenome: $o/panaroo/gene_presence_absence.csv)" >&2
