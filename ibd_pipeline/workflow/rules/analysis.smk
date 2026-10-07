# ---- analyses run for every (variant, layer): variant in {raw, mmuphin}, layer in {taxa, pathways}

def meta_of(wc):
    return f"{OUT}/merged/{wc.layer}/meta.tsv"

rule alpha_diversity:
    input: matrix=OUT + "/{variant}/taxa/matrix.tsv", meta=OUT + "/merged/taxa/meta.tsv"
    output:
        table=OUT + "/{variant}/taxa/alpha_diversity.csv",
        plot=OUT + "/{variant}/taxa/alpha_diversity.png",
    log: OUT + "/logs/alpha_{variant}.log"
    script: "../scripts/alpha_diversity.py"

rule beta_diversity:
    """Bray-Curtis, PCoA, PERMANOVA (batch + condition) and dispersion test."""
    input: matrix=OUT + "/{variant}/{layer}/matrix.tsv", meta=meta_of
    output:
        coords=OUT + "/{variant}/{layer}/pcoa_coords.csv",
        permanova=OUT + "/{variant}/{layer}/permanova.csv",
        dispersion=OUT + "/{variant}/{layer}/dispersion.csv",
    threads: 4
    log: OUT + "/logs/beta_{variant}_{layer}.log"
    script: "../scripts/beta_diversity.R"

rule beta_plot:
    input: coords=rules.beta_diversity.output.coords, permanova=rules.beta_diversity.output.permanova
    output: OUT + "/{variant}/{layer}/beta_diversity.png"
    log: OUT + "/logs/beta_plot_{variant}_{layer}.log"
    script: "../scripts/beta_plot.py"

rule differential_abundance:
    """Per-feature linear model on log10 abundance ~ condition (+ batch), BH-FDR."""
    input: matrix=OUT + "/{variant}/{layer}/matrix.tsv", meta=meta_of
    output:
        table=OUT + "/{variant}/{layer}/differential_abundance.csv",
        plot=OUT + "/{variant}/{layer}/differential_abundance.png",
    log: OUT + "/logs/da_{variant}_{layer}.log"
    script: "../scripts/differential_abundance.py"

rule taxonomic_set_enrichment:
    """fgsea of taxonomic groups (Phylum..Genus) ranked by the DA t-statistic."""
    input: OUT + "/{variant}/taxa/differential_abundance.csv"
    output:
        all=OUT + "/{variant}/taxa/sea_all_results.csv",
        sig=OUT + "/{variant}/taxa/sea_results.csv",
    log: OUT + "/logs/sea_{variant}.log"
    script: "../scripts/taxonomic_sea.R"

rule sea_plot:
    input: OUT + "/{variant}/taxa/sea_all_results.csv"
    output: OUT + "/{variant}/taxa/sea_plot.png"
    log: OUT + "/logs/sea_plot_{variant}.log"
    script: "../scripts/sea_plot.py"

rule focus_genera:
    input: matrix=OUT + "/{variant}/taxa/matrix.tsv", meta=OUT + "/merged/taxa/meta.tsv"
    output:
        stats=OUT + "/{variant}/taxa/focus_genera_stats.csv",
        plot=OUT + "/{variant}/taxa/focus_genera.png",
    log: OUT + "/logs/focus_{variant}.log"
    script: "../scripts/focus_genera.py"

rule compare_variants:
    """Raw vs MMUPHin: variance explained by batch/condition, PCoA, DA concordance."""
    input:
        raw_perm=OUT + "/raw/{layer}/permanova.csv", cor_perm=OUT + "/mmuphin/{layer}/permanova.csv",
        raw_pcoa=OUT + "/raw/{layer}/pcoa_coords.csv", cor_pcoa=OUT + "/mmuphin/{layer}/pcoa_coords.csv",
        raw_da=OUT + "/raw/{layer}/differential_abundance.csv", cor_da=OUT + "/mmuphin/{layer}/differential_abundance.csv",
    output:
        r2=OUT + "/comparison/{layer}/variance_explained.png",
        pcoa=OUT + "/comparison/{layer}/pcoa_before_after.png",
        da=OUT + "/comparison/{layer}/da_concordance.png",
        summary=OUT + "/comparison/{layer}/summary.csv",
    log: OUT + "/logs/compare_{layer}.log"
    script: "../scripts/compare_variants.py"
