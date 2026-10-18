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
| `scripts/01_download_genomes.sh` | reproduces the download |
| `scripts/plot_ctxm_plasmids_v2.py`, `config.toml`, `environment.yml` | the script under test, its config and conda env |
| `ctxm_plasmid_test.dochist.json` | dochist session log |

## Reproduce the download
```bash
atb fetch --tables assembly.parquet,run.parquet   # small tables; plain `atb fetch` builds a ~35 GB AMR index
bash scripts/01_download_genomes.sh
```

## Still to do before running the plot script
The script reads MOB-recon and AMRFinderPlus output, not raw assemblies, so these remain:
1. `mob_recon -i genomes/<SAMN>.fa.gz -o mobrecon/<SAMN>` per sample (needs `plasmid_<CLUSTER>.fasta` files and `mge.report.txt`).
2. `amrfinder -n <fasta> -o amrfinder/<SAMN>.amrfinder.tsv` per sample (needs the `Contig id` column to match the FASTA headers).
3. Optionally `metadata.tsv` (columns `sample`, `host`).
4. `conda env create -f environment.yml && python scripts/plot_ctxm_plasmids_v2.py --config config.toml`
