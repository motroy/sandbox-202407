# Test data for `plot_ctxm_plasmids_v2.py`

Short-read *E. coli* assemblies for testing the reference-guided blaCTX-M plasmid alignment script.

- **Source study:** [Antimicrob Agents Chemother 2018, AAC.01689-18](https://journals.asm.org/doi/full/10.1128/aac.01689-18); SRA table `input/SraRunInfo_6.csv` is BioProject **PRJNA379782** (266 *E. coli* BioSamples, Illumina HiSeq 2000, one run each).
- **Assemblies:** AllTheBacteria (ATB) short-read assemblies, fetched with [atb-cli](https://github.com/AllTheBacteria/atb-cli) v0.22.1 by BioSample. All 266 had an assembly, and all 266 downloaded (411 MB, `genomes/<SAMN>.fa.gz`; git-ignored). ATB assemblies are not the paper's own, so contig content may differ from the published ones.

## Layout
| Path | Content |
|---|---|
| `input/SraRunInfo_6.csv` | SRA run table (BioSample list) |
| `input/biosamples.txt` | 266 unique BioSample accessions |
| `input/atb_query.tsv` | ATB lookup: `sample_accession`, `aws_url` |
| `scripts/00_select_subset.py` | reproduces the subset selection from ATB's precomputed tables (needs `atb fetch --tables assembly.parquet,assembly_stats.parquet,mlst.parquet,amrfinderplus.parquet`, plus `duckdb`) |
| `scripts/01_download_genomes.sh` | reproduces the download |
| `scripts/02_make_subset.sh`, `scripts/03_run_mobrecon.sh` | link the subset; run MOB-recon (`mobrecon/`, git-ignored) |
| `scripts/plot_ctxm_plasmids_v2.py`, `config.toml`, `environment.yml` | the script under test, its config and conda env |
| `ctxm-plasmid-test.dochist.json` | dochist session log |

## Reproduce the download
```bash
atb fetch --tables assembly.parquet,run.parquet   # small tables; plain `atb fetch` builds a ~35 GB AMR index
bash scripts/01_download_genomes.sh
```

## Test subset (34 of 266 assemblies)
MOB-recon and plasmid clustering on all 266 genomes is unnecessary for testing, so `input/subset.tsv` lists 34 chosen for this script's use case (the script needs >=2 blaCTX-M isolates sharing a MOB-suite cluster, plus an optional non-CTX-M isolate to serve as the reference). Selection used ATB's precomputed AMRFinderPlus calls (`input/atb_amrfinderplus_prjna379782.tsv`), MLST and assembly stats:
- **27 blaCTX-M isolates** (65 of 266 carry it: CTX-M-9 x43, -15 x18, -27 x3, unspecified x2, -24 x1). For every variant + ST group the 2 best-assembled (fewest contigs) isolates are taken, so each probable shared backbone has >=2 members: CTX-M-9 in ST216 (x5, the largest group, spread over assembly quality), ST607, ST409, ST401, ST13850; CTX-M-15 in ST410, ST405, ST2363; CTX-M-27 in ST131. Singletons are kept only for the rare variants (CTX-M-24 ST38, CTX-M-27 ST3541, an unspecified blaCTX-M, ST2363) and the two dual-variant isolates (SAMN06630209, SAMN06630194) for the multi-gene case.
- **7 non-CTX-M controls**, the best-assembled isolate in each of the ST131/410/607/216/38/2363/401 groups, as candidate references.
- Clusters are **not known yet**: ST is only a proxy for a shared plasmid, so expect some of the intended MOB-suite clusters (AA147 etc.) to be skipped.
- `scripts/02_make_subset.sh` links the subset into `subset_genomes/` and copies ATB's AMRFinderPlus tables (contig IDs match the FASTA headers) into `amrfinder/`, which can replace step 2 below.

## Still to do before running the plot script
The script reads MOB-recon and AMRFinderPlus output, not raw assemblies, so these remain:
1. `mob_recon -i subset_genomes/<SAMN>.fa.gz -o mobrecon/<SAMN>` per sample (needs `plasmid_<CLUSTER>.fasta` files and `mge.report.txt`).
2. (optional for the subset, already staged by `02_make_subset.sh`) `amrfinder -n <fasta> -o amrfinder/<SAMN>.amrfinder.tsv` per sample (needs the `Contig id` column to match the FASTA headers).
3. Optionally `metadata.tsv` (columns `sample`, `host`).
4. `conda env create -f environment.yml && python scripts/plot_ctxm_plasmids_v2.py --config config.toml`
