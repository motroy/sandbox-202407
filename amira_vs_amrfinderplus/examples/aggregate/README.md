# Example: `amira_aggregate.py`

Re-runnable with only `python3`: `./run_example.sh` (regenerates `output/`).

**Inputs** (`input/`, copied from this repo's results)
- `ST131/amira_results.tsv`: Amira on SRR13116389 (100k ONT reads, *E. coli* ST131)
- `FWSEC0007/amira_results.tsv`: Amira on SRR7896269 (*E. coli* O103:H2, SAMN08768108)
- `*.amrfinderplus.tsv`: AMRFinderPlus on a same-isolate assembly (ST131: Flye on the same reads; FWSEC0007: dragonflye + Polypolish + pypolca)

**Outputs** (`output/`)
| file | what it is |
|---|---|
| `summary.tsv` | one row per sample: allele count, gene count, `est_copies`, mean identity/coverage/depth, gene list |
| `all_alleles.tsv` | every Amira allele row plus `gene_family` and `ambiguous_family` |
| `matrix_alleles.tsv` / `matrix_copies.tsv` | sample × gene family, allele counts / summed approximate copy number |
| `compare_amrfinderplus.tsv` | per gene family: `both`, `amira_only` or `amrfinder_only` with Amira alleles/copies and AMRFinderPlus hits |
| `summary_min_depth_0.7.tsv` | same summary after dropping alleles with relative depth < 0.7 (ST131: 72 → 48 alleles) |

Read allele counts with care: for ST131 they are inflated by read errors (e.g. blaTEM: 9 alleles, ~6.5 estimated copies). See the main README.
