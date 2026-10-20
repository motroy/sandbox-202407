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
`--style` adjusts the interactive plot through Playwright first: groups are coloured by class (blaCTX-M red, ISEcp1 blue, other transposases purple, other ARGs orange, rest grey), legend entries are renamed from the member genes' products instead of "Group N", only blaCTX-M / ISEcp1 / IS26 / ARG genes are labelled (with a white halo and extra row spacing), and a `<variant>.svg` and `<variant>.styled.html` are written next to the PNG. The gene classification is regex-based on the GBFF product text, so check it against your annotation.

The `flank_bp = 5000` default and the recurrence-based selection of representatives are my assumptions; the Methods don't give them. Clinker only aligns CDS with translations, so use annotated GBFFs.

Genomes sharing the same flanking-gene pair belong to the same integration region (group `insertion_sites.tsv` by `up_anchor_min`/`down_anchor_min`).

## Assembly quality
Use long-read or hybrid assemblies (the paper's gate: N50 >= 100 kb, completeness >= 90%; see `00_select_genomes.py filter`). IS elements and blaCTX-M are repeats, so short-read drafts break at them and chromosome/plasmid, flank and insertion-site calls are unreliable. `context_hits.tsv` has a `context_quality` column: `ok`, or `short_contig` (linear contig < `min_contig_bp`) and/or `truncated_flank` (linear contig ends within `window_bp` of the region). Circular contigs (MOB-recon `circularity_status`) are never flagged. Treat flagged hits as unreliable and consider excluding them from rates.

## Tools (chosen by the user; the paper used ABRicate/ResFinder/ISfinder/Platon instead)
| Task | Tool |
|---|---|
| ARGs, blaCTX-M variants | AMRFinderPlus (`--plus --organism Escherichia`; point mutations are excluded from barcodes) |
| Chromosome vs plasmid, plasmid clusters, IS/MGEs | MOB-recon (`contig_report.txt`, `mge.report.txt`) |
| ST / phylogroup / O:H / fimH | mlst, ClermonTyping, ECTyper, FimTyper |
| cgMLST de-duplication | chewBBACA AlleleCall + cgmlst-dists |
`05_typing.py` assigns ST131 subclades with a simple fimH / blaCTX-M-15 heuristic (fimH41 = A, fimH22 = B, fimH30 = C2 with blaCTX-M-15 else C1, other fimH = Other), not the paper's definitions (its ref. 47, which I have not seen). On the 30-genome test subset it matches the paper's clade for 11/11 ST131 genomes. ClermonTyping and FimTyper output layouts vary by version, so check `typing.tsv` after upgrades.

**Installing FimTyper** (not on conda): `git clone https://bitbucket.org/genomicepidemiology/fimtyper.git` and `.../fimtyper_db.git`, then `export FIMTYPER_DIR=<fimtyper clone> FIMTYPER_DB=<fimtyper_db clone>`. It is a Perl script needing BLAST, `perl-bioperl` and `perl-try-tiny-retry` (in `environment.yml`). ClermonTyping is the bioconda `clermontyping` package; its R report step may print an error, but `*_phylogroups.txt` is still written. Because of the MOB-recon choice the paper's Platon + 500 kb size rule and manual review of >500 kb contigs are not applied; gene-level and genome-level counts are in `summary_by_st.tsv`.

## Status
Re-run needed: the stored AMRFinderPlus outputs in `test_data/run/` were made without `--organism`; `run_subset.sh` skips existing files, so delete `test_data/run/amrfinder/` first. `context_quality` was checked on the existing outputs only.

`06` and `07` were run against the real `clinker` on two synthetic windows (one reverse-strand). `00`, `02`, `03` and `05` were run on small synthetic inputs only. `04_plot_sites.R` is **untested** (R was not available) and `01_run_tools.sh` has not been run.

## Not included
The manual SnapGene / ISfinder / TnCentral boundary curation of the most frequent insertions; the classification of shared IS into three transposition-mechanism categories (Supplementary Note 3, not available to me).
Genome screening/filtering (the paper started from 451,174 public genomes, keeping 2,269 CTX-M long-read/hybrid assemblies), ST/subclade typing, and the integration-rate modelling. Use short-read assemblies with care: contigs break at IS elements, so chromosome/plasmid calls are unreliable.

## Test data (Zenodo 10.5281/zenodo.20934804)
`test_data/summary_df_2507genomes.txt` is the paper's metadata table (2,507 genomes, with the authors' chromosomal/plasmidic blaCTX-M calls, ARG barcodes and ST131 clades, so it doubles as ground truth). Fetch it with `zenodo_get 20934804`.
```bash
python test_data/make_subset.py test_data/summary_df_2507genomes.txt test_data/subset.tsv   # 30 stratified genomes (seed 1)
python test_data/download_genomes.py test_data/subset.tsv test_data/genomes                  # FASTA + GBFF from NCBI Datasets
```
