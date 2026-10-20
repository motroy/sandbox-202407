# Example results (30-genome Zenodo test subset)

Final tables and Clinker plots only; raw AMRFinderPlus / MOB-recon / mlst / ECTyper output (about 150 MB) is not committed. Regenerate with `test_data/run_subset.sh`.

**Provenance:** AMRFinderPlus ran with `--plus` but *without* `--organism Escherichia` (added afterwards); the tables were regenerated from those outputs with the current scripts. ClermonTyping and FimTyper were run, so phylogroup, fimH and `st131_subclade` are filled (11/11 ST131 clades match the paper; the rule is a fimH/blaCTX-M-15 heuristic, see `scripts/05_typing.py`).

| File | One row per | Content |
|---|---|---|
| `context_hits.tsv` | blaCTX-M copy | contig, coordinates, chromosome/plasmid, plasmid cluster, relative position, inserted region, `context_quality`, ARG barcode + class, mediator, IS/MGEs and co-located ARGs in the window |
| `barcodes.tsv` | distinct barcode x location | recurrence over genomes, variants, STs, class vs plasmid barcodes |
| `summary_by_st.tsv` | ST | gene-level and genome-level counts, % chromosomal, mediators, % co-integrated |
| `insertion_sites.tsv` | chromosomal copy | K-12 position (min, 0-100), method, anchor-chain lengths |
| `site_histogram.tsv` | 0.1-min bin | genomes per insertion position |
| `typing.tsv` | genome | ST, phylogroup, O:H, fimH, ST131 subclade |
| `comparison.tsv` | genome | our calls next to the paper's (Zenodo table) |
| `clinker/<variant>.html/.csv/.png/.svg` | variant | synteny plot (interactive HTML, styled PNG/SVG from `scripts/07_clinker_png.py --style`: coloured gene classes, renamed legend, key genes labelled) and pairwise protein identities of the blaCTX-M windows; only variants with >= 2 windows are drawn |
