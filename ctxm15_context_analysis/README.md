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
3. Find IS/MGEs within `window_bp` of the gene from MOB-recon's `mge.report.txt` (`contig_start`/`contig_end`/`mge_type`) and assign the mediator, ISEcp1 first, then IS26.
4. List other AMR genes in the window (co-integration).
5. **ARG barcode** = sorted set of ARGs in the region. Chromosomal barcodes are compared with plasmid barcodes: `identical_to_plasmid`, `contained_in_larger_plasmid` (gene stabilised within part of a plasmid MDR island) or `chromosome_specific`. `barcodes.tsv` gives recurrence by genomes, variants and STs.
6. **Insertion sites** (`03_insertion_sites.py`): up to 5 CDS per flank are mapped to E. coli K-12 on the 0-100 min scale (gene symbol / locus tag / product). Per flank the longest contiguous anchor chain (gap <= 1.0 min, >= 2 genes) is kept; both flanks with junction gap <= 5.0 min give the circular midpoint, otherwise the longer flank or the single valid flank, else `ND`. Positions are rounded to 0.1 min; `04_plot_sites.R` draws the circular map (optional variant and ST filters).

## Usage
```bash
conda env create -f environment.yml
scripts/01_run_tools.sh SAMPLE1 assemblies/SAMPLE1.fa   # per sample
# fill input/samples.tsv (see input/samples.example.tsv; `st` optional, `gbff` needed for step 03)
python scripts/02_context_analysis.py --config config.toml
python scripts/05_typing.py config.toml               # typing.tsv + ST131 subclades
python scripts/03_insertion_sites.py config.toml      # needs input/k12_MG1655.gbff (NC_000913.3)
Rscript scripts/04_plot_sites.R results/insertion_sites.tsv results/sites.pdf blaCTX-M-15 ST131
```
Outputs (`results/`, git-ignored): `context_hits.tsv`, `barcodes.tsv`, `summary_by_st.tsv`, `insertion_sites.tsv`, `site_histogram.tsv`, `typing.tsv`, `clinker/`. A committed example from the test subset, with a description of every column set, is in [`test_data/example_results/`](test_data/example_results/).

**Clinker synteny (`06_clinker_windows.py`):** per variant, cuts the inserted region +- `flank_bp` from each GBFF, orients it so blaCTX-M is forward, and writes `results/clinker/<variant>/<sample>.gbk` plus a manifest. Genomes are ordered by barcode recurrence (most frequent structure first, one window per genome, capped by `max_per_variant`). `--run` then calls `clinker` to produce `<variant>.html` and `<variant>.csv`.
```bash
python scripts/06_clinker_windows.py --config config.toml --run
```
**PNG export (`07_clinker_png.py`):** `python scripts/07_clinker_png.py results/clinker` renders each `<variant>.html` to a cropped `<variant>.png` with headless Chromium (`pip install playwright`; pass `--browser` if Playwright's own Chromium is not installed).
`--style` adjusts the interactive plot through Playwright first: groups and their ribbons get pastel colours by class (blaCTX-M red, ISEcp1 blue, other transposases purple, other ARGs orange, rest light grey; `--grey-links` keeps the identity-shaded grey ribbons and identity bar instead), legend entries are renamed from the member genes' products instead of "Group N", only blaCTX-M / ISEcp1 / IS26 / ARG genes are labelled (with a white halo and extra row spacing), and a `<variant>.svg` and `<variant>.styled.html` are written next to the PNG. The gene classification is regex-based on the GBFF product text, so check it against your annotation.

The `flank_bp = 5000` default and the recurrence-based selection of representatives are my assumptions; the Methods don't give them. Clinker only aligns CDS with translations, so use annotated GBFFs.

Genomes sharing the same flanking-gene pair belong to the same integration region (group `insertion_sites.tsv` by `up_anchor_min`/`down_anchor_min`).

## Why ResFinder + Platon are the defaults (30-genome test subset vs the paper's table)
| | ResFinder + Platon | AMRFinderPlus + MOB-recon |
|---|---|---|
| blaCTX-M chromosome/plasmid copy numbers match | 30/30 | 28/30 |
| co-located ARG set (barcode partners) matches the paper | 8/8 exact | gene names differ (allele naming, `catB3`), window-bound |
ResFinder uses the paper's gene names (`aac(6')-Ib-cr`, `aac(3)-IIa`, `mph(A)`), found the CTX-M-172/183 copies AMRFinderPlus missed, and Platon matched the paper on the one contig where MOB-recon disagreed. Notes: ResFinder lists every tied allele as a separate hit (a novel allele can tie 5 ways), so overlapping hits are collapsed to one locus and the tied alleles are kept in `allele_ties`; ResFinder gives no strand (step 06 reads it from the GBFF). The ARG barcode is now the ARGs inside the *inserted region* (everything within `window_bp`, extended by ARGs/IS chained on with gaps <= `region_gap_bp`, capped at `region_max_bp`), which recovered islands longer than the flank window. Only 8 of the 30 genomes have a paper barcode, so treat the 8/8 as encouraging, not conclusive.

## IS detection benchmark (`11_compare_is_tools.py`, 29 chromosomal blaCTX-M copies, IS within 10 kb)
| Route | ISEcp1 | IS26 | none |
|---|---|---|---|
| ABRicate + ISfinder (ISEc9 = ISEcp1) | 23 (79%) | 5 (17%) | 1 |
| MOB-recon `mge.report.txt` (measured with `is_tool = "mobrecon"`) | 20 (69%) | 1 (3%) | 8 |
| ISEScan 1.7.3 (family: IS1380 / IS6) | 11 (38%) | 7 (24%) | 11 |
| ISEScan + BLAST naming vs ISfinder | 11 (38%) | 6 (21%) | 12 |
The paper (2,269 genomes) reports ISEcp1 77.3% and IS26 18.9%. ISfinder is the reference *database* and the only route here that matches those proportions; it is a plain BLAST search, so it gives no element boundaries and needs care with duplicate entries (ISEc9 = ISEcp1). ISEScan is a better general IS *inventory* tool (de novo, element boundaries) but missed 13 of the 23 ISEcp1 copies that ISfinder finds as full-length 100% matches (ISEcp1's inverted repeats are imperfect), so it is not safe for calling these two mediators. MOB-recon's MGE set has no IS26 entry. ISEScan was run on 60 kb windows around each hit (whole-genome runs took 15+ min each; `test_data/run_isescan_windows.py`).

## Phylogeny (`08_phylogeny.sh`, `09_tree_annotate.py`)
```bash
THREADS=8 scripts/08_phylogeny.sh test_data/genomes results/phylo      # <id>.fa + <id>.gbff per genome
python scripts/09_tree_annotate.py config.toml results/phylo
```
1. **Mashtree** on the assemblies (`mashtree.dnd`).
2. **Panaroo** on the NCBI GBFFs (no re-annotation): strict clean mode, core threshold 0.98, `--remove-invalid-genes` (NCBI annotations contain frameshifted CDS that Panaroo rejects), FAMSA with `--codons` (protein alignment back-translated; about 5x faster than MAFFT here, `ALIGNER=mafft` to switch). Panaroo only converts GenBank files when `-i` is a text file listing the paths, which the script does.
3. **IQ-TREE 3** on the core-gene alignment with a fixed `GTR+F+I+G4` (`MODEL=MFP` runs ModelFinder, which takes hours on a multi-Mb alignment). Support: ultrafast bootstrap by default (`BOOTSTRAPS=1000`, about an hour for 30 genomes; `BOOTSTRAPS=0 ALRT=1000` uses fast SH-aLRT instead).
4. `09_tree_annotate.py` writes `annotation/`: `tree_metadata.tsv` (ST, phylogroup, subclade, variant, integration status, mediators), iTOL colour-strip files, `tree_<name>.png` and `tree_<name>_cladogram.png` (midpoint-rooted, strips for integration / ST / variant), and the normalised Robinson-Foulds distance between the two topologies.

The core threshold, clean mode, model and support settings are my choices; the paper used PPanGGOLiN and FastTree instead. Install: `panaroo`, `iqtree`, `mashtree`, `famsa` (and `mafft`) from bioconda (kept in a separate environment from the AMR tools if conda resolution complains); `pip install matplotlib biopython` for step 09.

## Comparison with the paper's own tools (`test_data/run_paper_tools.sh`, `scripts/10_compare_paper_tools.py`)
Runs the tools named in the paper's Methods on the same genomes and compares them with this pipeline and with the paper's table:
- **ABRicate 1.0.1** with `resfinder`, `ecoh` (bundled) and two custom databases you build once: `isfinder` (`IS.fna` from github.com/thanhleviet/ISfinder-sequences, the paper used the Oct 2020 ISfinder; this copy's date is unknown) and `fimtyper` (`fimH.fsa` from bitbucket.org/genomicepidemiology/fimtyper_db): `mkdir $CONDA_PREFIX/db/<name>`, copy the FASTA there as `sequences`, run `abricate --setupdb`.
- **Platon 1.7** (database from zenodo.org/records/4066768, 1.7 GB) with the paper's size rule (contigs > 500 kb are chromosomal; a Platon plasmid call on one is flagged). Platon writes no table when no contig passes its filters, which the comparison reads as "no plasmid contigs".
- **PPanGGOLiN 2.2.4** (`workflow` on the GBFFs, then `msa --partition softcore --soft_core 0.95 --source dna --phylo`) and **FastTree 2.1.10** (`-gtr -gamma -nt`) for the paper's phylogeny.
`10_compare_paper_tools.py <config> <subset.tsv> <run_paper_dir> [fasttree.nwk]` writes `comparison/`: `ctxm_hits_paper_tools.tsv`, `per_genome_comparison.tsv`, `tree_rf.tsv` (normalised Robinson-Foulds between FastTree, IQ-TREE and Mashtree) and `summary.txt`. Caveats: ISfinder lists ISEc9 as a 100%-identical copy of ISEcp1 and ABRicate prints only one name, so ISEc9 is counted as ISEcp1; MOB-suite's MGE set has no IS26 entry, so this pipeline counts its IS6-family elements as IS26-family.

## ggCaller vs NCBI annotations (30 genomes; `ggcaller` 1.5.2 gene calls, then Panaroo)
In 1.5.2 ggCaller does graph-based gene calling and ORF clustering only (no pangenome table), so the comparison is Panaroo on ggCaller's calls vs Panaroo on the NCBI annotations (same settings: strict, `--remove-invalid-genes`, FAMSA codons, core 0.98, then IQ-TREE `GTR+F+I+G4` + SH-aLRT).
| | NCBI annotations | ggCaller calls |
|---|---|---|
| core (>= 99%) / soft-core / shell / cloud | 2,930 / 437 / 2,825 / 4,965 | 3,039 / 252 / 2,824 / 5,012 |
| total gene families | 11,157 | 11,127 |
| parsimony-informative core sites | 122,655 | 144,340 |
Same gene calls for every genome convert soft-core families into core ones (fewer annotation-dependent gaps). The two IQ-TREE topologies differ by normalised RF 0.41 (each differs from Mashtree by 0.48-0.52). The 15 branches that are well supported (SH-aLRT >= 95) in both trees include the backbone between STs; the 12 well-supported branches that contradict the other tree all sit in near-identical groups (the ST131 cluster and three ST38 genomes). SH-aLRT is therefore overconfident there: for near-clonal groups use ultrafast bootstrap and mask recombination (Gubbins/ClonalFrameML) before reading fine structure. ggCaller needs its 120 MB model database (`ftp.ebi.ac.uk/pub/databases/pp_dbs/ggCallerdb.tar.bz2`) unpacked into the package's `models/` directory (the built-in downloader hung behind a proxy here). Whether ggCaller is the better default at 2,500 genomes is untested here.

## Assembly quality
Use long-read or hybrid assemblies (the paper's gate: N50 >= 100 kb, completeness >= 90%; see `00_select_genomes.py filter`). IS elements and blaCTX-M are repeats, so short-read drafts break at them and chromosome/plasmid, flank and insertion-site calls are unreliable. `context_hits.tsv` has a `context_quality` column: `ok`, or `short_contig` (linear contig < `min_contig_bp`) and/or `truncated_flank` (linear contig ends within `window_bp` of the region). Circular contigs (MOB-recon `circularity_status`) are never flagged. Treat flagged hits as unreliable and consider excluding them from rates.

## Tools (chosen by the user; the paper used ABRicate/ResFinder/ISfinder/Platon instead)
| Task | Tool |
|---|---|
| ARGs, blaCTX-M variants, barcodes | **ResFinder 4.7** (default; 80% identity/coverage) or AMRFinderPlus (`arg_tool = "amrfinder"`; `--plus --organism Escherichia`, point mutations excluded) |
| Chromosome vs plasmid | **Platon 1.7 + the paper's 500 kb size rule** (default) or MOB-recon (`location_tool = "mobrecon"`); if both are supplied, disagreements are flagged in `location_disagreement` |
| IS mediators (ISEcp1, IS26) | **ABRicate vs ISfinder** (`is_tool = "isfinder"`, default; ISEc9 is counted as ISEcp1) or MOB-recon MGEs (`is_tool = "mobrecon"`) |
| Plasmid clusters | MOB-recon `contig_report.txt` (optional) |
| ST / phylogroup / O:H / fimH | mlst, ClermonTyping, ECTyper, FimTyper |
| cgMLST de-duplication | chewBBACA AlleleCall + cgmlst-dists |
| Phylogeny | Mashtree (assemblies), Panaroo (pangenome + core alignment, FAMSA), IQ-TREE 3 |
`05_typing.py` assigns ST131 subclades with a simple fimH / blaCTX-M-15 heuristic (fimH41 = A, fimH22 = B, fimH30 = C2 with blaCTX-M-15 else C1, other fimH = Other), not the paper's definitions (its ref. 47, which I have not seen). On the 30-genome test subset it matches the paper's clade for 11/11 ST131 genomes. ClermonTyping and FimTyper output layouts vary by version, so check `typing.tsv` after upgrades.

**Installing FimTyper** (not on conda): `git clone https://bitbucket.org/genomicepidemiology/fimtyper.git` and `.../fimtyper_db.git`, then `export FIMTYPER_DIR=<fimtyper clone> FIMTYPER_DB=<fimtyper_db clone>`. It is a Perl script needing BLAST, `perl-bioperl` and `perl-try-tiny-retry` (in `environment.yml`). ClermonTyping is the bioconda `clermontyping` package; its R report step may print an error, but `*_phylogroups.txt` is still written. Because of the MOB-recon choice the paper's Platon + 500 kb size rule and manual review of >500 kb contigs are not applied; gene-level and genome-level counts are in `summary_by_st.tsv`.

## Status
Re-run needed: the stored AMRFinderPlus outputs in `test_data/run/` were made without `--organism`; `run_subset.sh` skips existing files, so delete `test_data/run/amrfinder/` first. `context_quality` was checked on the existing outputs only.

`08` and `09` were run on the 30-genome subset (2,930 core genes, 2.75 Mb core alignment, SH-aLRT 1000, no ultrafast bootstrap). `06` and `07` were run against the real `clinker` on two synthetic windows (one reverse-strand). `00`, `02`, `03` and `05` were run on small synthetic inputs only. `04_plot_sites.R` is **untested** (R was not available) and `01_run_tools.sh` has not been run.

## Not included
The manual SnapGene / ISfinder / TnCentral boundary curation of the most frequent insertions; the classification of shared IS into three transposition-mechanism categories (Supplementary Note 3, not available to me).
Genome screening/filtering (the paper started from 451,174 public genomes, keeping 2,269 CTX-M long-read/hybrid assemblies), ST/subclade typing, and the integration-rate modelling. Use short-read assemblies with care: contigs break at IS elements, so chromosome/plasmid calls are unreliable.

## Test data (Zenodo 10.5281/zenodo.20934804)
`test_data/summary_df_2507genomes.txt` is the paper's metadata table (2,507 genomes, with the authors' chromosomal/plasmidic blaCTX-M calls, ARG barcodes and ST131 clades, so it doubles as ground truth). Fetch it with `zenodo_get 20934804`.
```bash
python test_data/make_subset.py test_data/summary_df_2507genomes.txt test_data/subset.tsv   # 30 stratified genomes (seed 1)
python test_data/download_genomes.py test_data/subset.tsv test_data/genomes                  # FASTA + GBFF from NCBI Datasets
```
