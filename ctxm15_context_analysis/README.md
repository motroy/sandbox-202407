# blaCTX-M-15 genomic-context analysis

Modelled on Danjean et al., *Chromosomal integration of blaCTX-M redefines the evolutionary trajectory of antibiotic resistance in E. coli through lineage-structured expansion* ([bioRxiv 2026.10.07.757360](https://www.biorxiv.org/content/10.64898/2026.10.07.757360v1.full)).

**Caveat:** bioRxiv rate-limited the full text, so this is built from the abstract only. The paper's exact tools, thresholds and statistics (subclade modelling, hotspot analysis) are **not** reproduced; check the Methods and adjust `config.toml`.

## What it does (per blaCTX-M-15 hit, long-read/hybrid assemblies)
1. Locate the gene (AMRFinderPlus) and classify its contig as chromosome or plasmid (MOB-recon).
2. Record the relative position along the contig (the paper reports integration across the whole chromosome).
3. Find IS elements within `window_bp` of the gene (BLAST vs an IS reference set) and assign the mediator, ISEcp1 first, then IS26.
4. List other AMR genes in the window (co-integration).
5. Summarise by ST: % chromosomal, mediator counts, % co-integrated.

## Usage
```bash
conda env create -f environment.yml
scripts/01_run_tools.sh SAMPLE1 assemblies/SAMPLE1.fa is_references.fa   # per sample
# fill input/samples.tsv (see input/samples.example.tsv; `st` is optional)
python scripts/02_context_analysis.py --config config.toml
```
Outputs: `results/context_hits.tsv`, `results/summary_by_st.tsv`.

## Not included
Genome screening/filtering (the paper started from 451,174 public genomes, keeping 2,269 CTX-M long-read/hybrid assemblies), ST/subclade typing, and the integration-rate modelling. Use short-read assemblies with care: contigs break at IS elements, so chromosome/plasmid calls are unreliable.
