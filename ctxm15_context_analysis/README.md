# blaCTX-M-15 genomic-context analysis

Modelled on Danjean et al., *Chromosomal integration of blaCTX-M redefines the evolutionary trajectory of antibiotic resistance in E. coli through lineage-structured expansion* ([bioRxiv 2026.10.07.757360](https://www.biorxiv.org/content/10.64898/2026.10.07.757360v1.full)).

**Caveat:** built from the abstract plus the Methods excerpts for ARG barcodes, representative structures and insertion-site mapping. The earlier Methods (genome screening, annotation tools, thresholds, which are in the paper's Supplementary Information) are not available to me, so the window size, identity cutoffs and tool choices in `config.toml` are my assumptions.

## Step 0: genome collection (`00_select_genomes.py`, `00b_cgmlst.sh`)
1. `filter`: from an NCBI assembly metadata table keep taxid 562, ANI `OK`, N50 >= 100 kb, completeness >= 90%, contamination <= 5%, BioProject and release date present, long-read/hybrid technology (column names and thresholds in `[selection]`). Download FASTA + GBFF for these.
2. Run ResFinder per genome (`<accession>.tsv`), then `00b_cgmlst.sh` (chewBBACA AlleleCall, Achtman 2,513-gene scheme, cgmlst-dists).
3. `dedup`: keep genomes with a blaCTX-M hit (ResFinder identity and coverage >= 80%), then collapse duplicates (same BioProject AND cgMLST distance <= 10 AND identical ARG profile) to one representative per cluster. The representative rule (highest N50, then completeness) is my assumption, as the Methods don't say.
```bash
python scripts/00_select_genomes.py filter ncbi_meta.tsv candidates.tsv
python scripts/00_select_genomes.py dedup candidates.tsv resfinder/ cgmlst_dists.tsv final_genomes.tsv
```
The Methods text I have stops at "Phylogroups were assigned using ClermontTyper", so phylogroup / ST typing is not implemented. The downstream scripts use AMRFinderPlus for ARG calls, whereas the paper used ResFinder for selection; the paper's downstream ARG tool is not known to me.

## What steps 1-6 do (per blaCTX-M-15 hit, long-read/hybrid assemblies)
1. Locate the gene (AMRFinderPlus) and classify its contig as chromosome or plasmid (MOB-recon).
2. Record the relative position along the contig (the paper reports integration across the whole chromosome).
3. Find IS elements within `window_bp` of the gene (BLAST vs an IS reference set) and assign the mediator, ISEcp1 first, then IS26.
4. List other AMR genes in the window (co-integration).
5. **ARG barcode** = sorted set of ARGs in the region. Chromosomal barcodes are compared with plasmid barcodes: `identical_to_plasmid`, `contained_in_larger_plasmid` (gene stabilised within part of a plasmid MDR island) or `chromosome_specific`. `barcodes.tsv` gives recurrence by genomes, variants and STs.
6. **Insertion sites** (`03_insertion_sites.py`): up to 5 CDS per flank are mapped to E. coli K-12 on the 0-100 min scale (gene symbol / locus tag / product). Per flank the longest contiguous anchor chain (gap <= 1.0 min, >= 2 genes) is kept; both flanks with junction gap <= 5.0 min give the circular midpoint, otherwise the longer flank or the single valid flank, else `ND`. Positions are rounded to 0.1 min; `04_plot_sites.R` draws the circular map (optional variant and ST filters).

## Usage
```bash
conda env create -f environment.yml
scripts/01_run_tools.sh SAMPLE1 assemblies/SAMPLE1.fa is_references.fa   # per sample
# fill input/samples.tsv (see input/samples.example.tsv; `st` optional, `gbff` needed for step 03)
python scripts/02_context_analysis.py --config config.toml
python scripts/03_insertion_sites.py config.toml      # needs input/k12_MG1655.gbff (NC_000913.3)
Rscript scripts/04_plot_sites.R results/insertion_sites.tsv results/sites.pdf blaCTX-M-15 ST131
```
Outputs: `results/context_hits.tsv`, `barcodes.tsv`, `summary_by_st.tsv`, `insertion_sites.tsv`, `site_histogram.tsv`.

Genomes sharing the same flanking-gene pair belong to the same integration region (group `insertion_sites.tsv` by `up_anchor_min`/`down_anchor_min`).

## Status
`00`, `02` and `03` were run on small synthetic inputs only. `04_plot_sites.R` is **untested** (R was not available) and `01_run_tools.sh` has not been run.

## Not included
Clinker synteny of flanking windows per variant (`clinker windows/*.gbk -p`) and the manual SnapGene / ISfinder / TnCentral boundary curation of the most frequent insertions; the classification of shared IS into three transposition-mechanism categories (Supplementary Note 3, not available to me).
Genome screening/filtering (the paper started from 451,174 public genomes, keeping 2,269 CTX-M long-read/hybrid assemblies), ST/subclade typing, and the integration-rate modelling. Use short-read assemblies with care: contigs break at IS elements, so chromosome/plasmid calls are unreliable.
