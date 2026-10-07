# IBD case/control pipeline (Snakemake)

A configurable IBD analysis built from the `curatedmetagenomicdata_snakemake` notebook reproduction. It differs in four ways:

1. **IBD only**: Healthy vs IBD, taxonomic and functional layers.
2. **Dataset selection**: pick any studies from the `curatedMetagenomicData` package, and/or add your own **local MetaPhlAn4 and HUMAnN outputs**. Both can be combined in one run.
3. **Batch-effect correction is optional**: run **without** it (`raw`), **with** MMUPHin (`mmuphin`), or **both** with an automatic comparison.
4. The condition is defined by explicit case/control rules (see [Case/control definition](#casecontrol-definition)).

## Quick start

```bash
cd ibd_pipeline
pixi install
pixi run setup                          # one-off fix-ups, incl. compiling MMUPHin (~20-30 min first time)
pixi run run                            # default config: 5 curated IBD studies, raw + MMUPHin
# curated + local data (uses the test fixtures in tests/fixtures):
pixi run snakemake --profile profiles/default --configfile config/config.example_local.yaml
```

`--configfile` overrides only the keys it contains; everything else comes from `config/config.yaml`. You can also override single values, e.g. `--config batch_correction='{mode: none}'`.

## Choosing datasets

All options are in `config/config.yaml`.

### From curatedMetagenomicData

```yaml
datasets:
  curated:
    studies: [HMP_2019_ibdmdb, HallAB_2017, NielsenHB_2014, LiJ_2014, IjazUZ_2017]   # any study names
    body_site: stool
    condition_column: study_condition
    control_regex: "^(control|healthy)$"
    case_regex: "^(IBD|UC|CD|ulcerative_colitis|crohns_disease)$"
```

Each study is downloaded once (latest snapshot) as `relative_abundance` (MetaPhlAn species) and `pathway_abundance` (HUMAnN pathways). Files are cached in `resources/experimenthub_cache/`.

### From local MetaPhlAn4 / HUMAnN output

```yaml
datasets:
  local:
    - name: my_cohort                      # unique dataset name; used as the default batch label
      metaphlan: "path/to/metaphlan/*_profile.txt"             # glob of per-sample profiles, OR one merged table
      humann_pathabundance: "path/to/humann/*_pathabundance.tsv"   # optional; glob OR merged (humann_join_tables)
      metadata: path/to/metadata.tsv       # tab-separated (.csv = comma-separated)
      sample_id_column: sample_id
      condition_column: diagnosis
      control_values: [healthy]            # exact (case-insensitive) values ...
      case_values: [CD, UC]
      # ... or regexes: control_regex / case_regex
      # sample_id_regex: "_profile$"       # optional: regex stripped from file names / column headers
```

Accepted formats:

| Input | Format |
|---|---|
| `metaphlan` | MetaPhlAn 3/4 per-sample profiles (`#clade_name  NCBI_tax_id  relative_abundance ...`), and/or tables from `merge_metaphlan_tables.py`. Only species-level (`s__`) clades are used; SGB (`t__`) lines and higher ranks are ignored. |
| `humann_pathabundance` | HUMAnN `*_pathabundance.tsv` (per sample) or merged tables. Stratified rows (`PWY|g__...`) are dropped, as are `UNMAPPED` / `UNINTEGRATED`. |
| `metadata` | One row per sample with the sample id, condition and, optionally, any other columns. |

Sample ids are taken from the file names (suffixes like `_profile.txt` and `_pathabundance.tsv` are stripped), or from the column headers of merged tables, and must match `sample_id_column`. Samples without metadata are dropped with a warning. If nothing matches, the error shows both id styles so that you can set `sample_id_regex`.

### Case/control definition

Samples are labelled `Healthy` (control) or `IBD` (case) (`condition_labels`) from your regexes/values; **samples matching neither are dropped**. For example, `LiJ_2014` also contains type-1 and type-2 diabetes samples; they are excluded by the default `case_regex`.

### How datasets are combined

- **Taxa** are matched on the species name (`s__...`), not the full lineage, so different MetaPhlAn database versions agree on most species. The lineage of the first dataset that contains a species is kept for higher-rank analyses. If you mix **MetaPhlAn 4 (SGB-era) and curatedMetagenomicData's MetaPhlAn 3 profiles**, species naming differs for some taxa, so check the number of shared features in `results/logs/merge_datasets.log`. `merge.join: inner` keeps shared species only; `outer` keeps the union and fills 0, which can bias comparisons across databases.
- **Pathways** (HUMAnN) are matched on pathway name. Datasets without pathway data are skipped for that layer, and `results/merged/sample_counts.csv` lists the samples that took part.
- Every sample is rescaled to sum to 100 after joining (`merge.renormalize`).

## Batch correction

```yaml
batch_correction:
  mode: both            # none | mmuphin | both
  batch_column: dataset # any metadata column present for all samples
```

- `none`: analyses on the merged data only (`results/raw/`).
- `mmuphin`: analyses on MMUPHin `adjust_batch`-corrected data (`results/mmuphin/`), with the batch as the batch variable and **condition as a protected covariate**.
- `both`: both of the above plus `results/comparison/`, showing how much variance the batch explained before and after, a PCoA before/after, and the concordance of the differential-abundance statistics.

MMUPHin needs at least two batches, each with at least `min_batch_size` samples. Features present in under `min_prevalence` of the samples are dropped before correction. Running on a single dataset therefore requires `mode: none`.

## Outputs (`results/`)

| Path | Contents |
|---|---|
| `merged/sample_counts.csv` | samples per dataset, condition and layer |
| `{raw,mmuphin}/{taxa,pathways}/matrix.tsv` | analysed abundance matrix (percent) |
| `.../beta_diversity.png`, `permanova.csv`, `dispersion.csv`, `pcoa_coords.csv` | Bray-Curtis PCoA, PERMANOVA (`~ batch + condition`, marginal; and `~ condition`), beta dispersion |
| `.../differential_abundance.{csv,png}` | linear model on log10 abundance ~ condition (optionally + batch), BH-FDR |
| `{raw,mmuphin}/taxa/alpha_diversity.{csv,png}` | Shannon, Simpson, richness, with Mann-Whitney tests |
| `{raw,mmuphin}/taxa/sea_*.csv`, `sea_plot.png` | fgsea enrichment of Phylum..Genus groups |
| `{raw,mmuphin}/taxa/focus_genera*` | focused genus/species plots (`analysis.focus_genera`) |
| `comparison/{taxa,pathways}/` | raw vs MMUPHin: `variance_explained.png`, `pcoa_before_after.png`, `da_concordance.png`, `summary.csv` |
| `logs/` | one log per job |

## Differences from the notebook pipeline

- **No Faith's PD.** It needs a phylogenetic tree, which local MetaPhlAn4 data don't come with.
- **HUMAnN pathways** are analysed: curatedMetagenomicData's data type is `pathway_abundance` (the notebook queried a non-existent `pathabundance` type, so it always reported "no pathways").
- **Case/control rule**: the notebook labelled every non-control sample "IBD", which put the 110 T1D/T2D samples of `LiJ_2014` in the IBD group. Here the case definition is explicit.
- The notebook's PERMANOVA used only `condition`. Here the batch is included as a second term, so the batch and condition effects can be compared.

## Layout

```
Snakefile  config/  profiles/  setup_env.sh  pixi.toml
workflow/rules/{data,correction,analysis}.smk
workflow/scripts/   R (fetch, MMUPHin, vegan, fgsea) and Python (parsing, stats, plots)
tests/              make_fixtures.R + fixtures (real curatedMetagenomicData samples re-exported as local MetaPhlAn/HUMAnN files)
```

See `setup_env.sh` for the environment fix-ups (rbiom 1.0.3 pin, curatedMetagenomicData payload, MMUPHin from Bioconductor).
