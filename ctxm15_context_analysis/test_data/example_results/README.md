# Example results (30-genome Zenodo test subset)

Final tables, plots and trees only; raw tool output (about 8 GB) is not committed. Regenerate with `test_data/run_subset.sh` (see the main README for the tool databases).

**Provenance.** Step 02 ran with the default callers: ResFinder 4.7 (80% identity/coverage), Platon 1.7 + the 500 kb size rule, ABRicate/ISfinder for IS elements (ISEc9 counted as ISEcp1). Typing: mlst, ClermonTyping, ECTyper, FimTyper. Plasmid/chromosome and IS calls are from this run; the older AMRFinderPlus/MOB-recon results are not included (the main README reports the comparison). ST131 subclades use the fimH / blaCTX-M-15 heuristic of `scripts/05_typing.py`.

| File | One row per | Content |
|---|---|---|
| `context_hits.tsv` | blaCTX-M copy | contig, coordinates, chromosome/plasmid, `location_disagreement`, `allele_ties` (ResFinder alleles tied for the locus), inserted region, `context_quality`, ARG barcode + class, mediator, IS and co-located ARGs |
| `barcodes.tsv` | distinct barcode x location | recurrence over genomes, variants, STs, class vs plasmid barcodes |
| `summary_by_st.tsv` | ST | gene-level and genome-level counts, % chromosomal, mediators, % co-integrated |
| `insertion_sites.tsv`, `site_histogram.tsv` | chromosomal copy; 0.1-min bin | K-12 position (min, 0-100), method, anchor chains |
| `typing.tsv` | genome | ST, phylogroup, O:H, fimH, ST131 subclade |
| `comparison.tsv` | genome | our calls next to the paper's (Zenodo table), incl. barcode partner Jaccard |
| `clinker/<variant>.html/.csv/.png/.svg` | variant | synteny of the blaCTX-M windows (`07_clinker_png.py --style`: pastel gene/ribbon classes, renamed legend, key genes labelled); only variants with >= 2 windows are drawn |
| `phylo/` | | `mashtree.dnd`; `iqtree.treefile` (Panaroo on NCBI annotations, GTR+F+I+G4, SH-aLRT 1000); `iqtree_ggcaller.treefile` (Panaroo on ggCaller gene calls, same settings); the Panaroo summaries; `annotation/` = tree PNGs (incl. cladograms), iTOL strips, `tree_metadata.tsv` |
| `paper_tools/` | | the paper's own tools (ABRicate ResFinder/ISfinder/FimTyper/EcOH, Platon) vs this pipeline (`10_compare_paper_tools.py`) and the IS-route benchmark (`11_compare_is_tools.py`) |

Headline numbers (30 genomes): copy numbers match the paper table 30/30, ST 30/30, ST131 clades 11/11, barcode partners 8/8 exact (only 8 genomes carry a paper barcode).
