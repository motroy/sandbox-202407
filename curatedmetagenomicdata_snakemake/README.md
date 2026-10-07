# curatedMetagenomicData – Snakemake pipeline

Reproduces the Colab notebook *"Using curatedMetagenomicData in colab"* as a
Snakemake workflow. The notebook mixed Python and R through `rpy2`, with state shared
between cells; here each analysis step is a standalone R or Python script, and the
hand-offs are explicit files.

## Run

```bash
cd curatedmetagenomicdata_snakemake
curl -fsSL https://pixi.sh/install.sh | bash   # if pixi is not installed (see ../install.pixi.sh)
pixi install                                    # R 4.5 + Bioconductor + Python + Snakemake (conda-forge/bioconda)
pixi run setup                                  # one-off fix-ups, see setup_env.sh
pixi run run                                    # = snakemake --profile profiles/default
```

Everything lands in `results/` (logs in `results/logs/`). Downloaded datasets are cached in
`resources/experimenthub_cache/`. Parameters (studies, thresholds, seeds) are in `config/config.yaml`.

## Workflow ↔ notebook

| Section | Rules | Notebook cells | Main outputs (`results/…`) |
|---|---|---|---|
| Demo (AsnicarF_2017) | `demo_asnicar`, `demo_plot` | 4–7 | `demo/demo_pca_alpha.png`, `demo/demo_stats.txt` (PERMANOVA, lm) |
| Alcohol & stool microbiota | `alcohol_study`, `alcohol_seaborn` | 9–11 | `alcohol/alcohol_ggplot_*.png`, `alcohol/alcohol_seaborn_summary.png` |
| IBD data assembly | `ibd_prepare` | 13–16 | `ibd/data_availability.csv`, `ibd/ibd_tse.rds`, `ibd/ibd_meta_alpha.csv`, `ibd/ibd_pcoa_coords.csv`, `ibd/ibd_bray_dist.rds` |
| IBD alpha/beta | `ibd_alpha_beta_plots`, `ibd_dist_to_csv`, `ibd_beta_stats` | 17–19 | `ibd/ibd_alpha_diversity.png`, `ibd/ibd_alpha_beta_overview.png`, `ibd/ibd_beta_statistical_summary.png`, `ibd/ibd_permanova.txt` |
| Taxonomic set enrichment | `ibd_taxonomic_sea`, `ibd_taxonomic_sea_plot` | 20–21 | `ibd/ibd_taxonomic_sea_results.csv`, `…_plot.png` |
| Marker clade enrichment | `ibd_marker_sea`, `ibd_marker_sea_plot` | 22–24 | `ibd/ibd_marker_clade_sea_results.csv`, `…_plot.png` |
| Coprococcus | `ibd_coprococcus_extract`, `ibd_coprococcus_stats` | 25–26 | `ibd/ibd_coprococcus_species_stats.csv`, `…_summary_plot.png` |
| MMUPHin batch correction (new, not in notebook) | `ibd_mmuphin`, `ibd_mmuphin_plot` | – | `ibd/mmuphin_pcoa_before_after.png`, `ibd/mmuphin_variance_explained.png`, `ibd/mmuphin_permanova.csv`, `ibd/mmuphin_adjusted_abundance.csv` |

Cells 0–2 (apt/pip installs, Google Drive mounting) are replaced by the pixi environment.

## Deviations from the notebook

- **Environment**: pixi (conda-forge + bioconda) instead of apt/pip/BiocManager on Colab. Versions come from the
  solver (R 4.5, Bioconductor 3.22: mia 1.18, curatedMetagenomicData 3.18.0).
- **`setup_env.sh`**: `curatedMetagenomicData` needs `rbiom::unifrac`, which was removed from current rbiom, so rbiom 1.0.3 is installed
  from the CRAN archive. The bioconda data package also installs its payload in a post-link script that pixi skips, so it is run manually.
- **Demo PCA**: `runPCA()` needs a `logcounts` assay, which the data lacks, so `logNormCounts()` is run first.
- **Serialised downloads**: ExperimentHub's cache index is not safe for concurrent writers, so jobs that download data share
  one `ehub` resource (set in `profiles/default/config.yaml`).
- **PERMANOVA**: the notebook built the skbio distance matrix from all samples but the grouping table from the
  NA-filtered subset; here the matrix is restricted to the same samples.
- **Outputs**: full (unfiltered) enrichment tables are also saved (`*_all_results.csv`); the files with the notebook's names hold
  only FDR < 0.05 rows. Plots are written as PNG files rather than displayed inline.
- Random-dependent steps (UMAP, fgsea, permutation tests) are seeded from `config/config.yaml`.
