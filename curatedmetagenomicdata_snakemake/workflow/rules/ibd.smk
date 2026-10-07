IBD = f"{OUT}/ibd"

rule ibd_prepare:
    """Load + merge the 5 IBD studies; alpha diversity, Bray-Curtis, PCoA."""
    output:
        availability=f"{IBD}/data_availability.csv",
        tse=f"{IBD}/ibd_tse.rds",
        meta=f"{IBD}/ibd_meta_alpha.csv",
        dist=f"{IBD}/ibd_bray_dist.rds",
        pcoa=f"{IBD}/ibd_pcoa_coords.csv",
    resources: ehub=1  # ExperimentHub cache is not safe for concurrent writers
    log: f"{OUT}/logs/ibd_prepare.log"
    script: "../scripts/ibd_prepare.R"

rule ibd_alpha_beta_plots:
    input:
        meta=rules.ibd_prepare.output.meta,
        pcoa=rules.ibd_prepare.output.pcoa,
        dist=rules.ibd_prepare.output.dist,
    output:
        overview=f"{IBD}/ibd_alpha_beta_overview.png",
        alpha=f"{IBD}/ibd_alpha_diversity.png",
        plot_df=f"{IBD}/ibd_plot_samples.csv",
    log: f"{OUT}/logs/ibd_alpha_beta_plots.log"
    script: "../scripts/ibd_alpha_beta_plots.py"

rule ibd_beta_stats:
    input:
        plot_df=rules.ibd_alpha_beta_plots.output.plot_df,
        pcoa=rules.ibd_prepare.output.pcoa,
        dist=f"{IBD}/ibd_bray_dist.csv",
    output:
        fig=f"{IBD}/ibd_beta_statistical_summary.png",
        txt=f"{IBD}/ibd_permanova.txt",
    log: f"{OUT}/logs/ibd_beta_stats.log"
    script: "../scripts/ibd_beta_stats.py"

rule ibd_dist_to_csv:
    """Export the Bray-Curtis distance matrix (RDS) as CSV for the Python steps."""
    input: rules.ibd_prepare.output.dist
    output: f"{IBD}/ibd_bray_dist.csv"
    log: f"{OUT}/logs/ibd_dist_to_csv.log"
    script: "../scripts/dist_to_csv.R"

rule ibd_taxonomic_sea:
    input: rules.ibd_prepare.output.tse
    output:
        all=f"{IBD}/ibd_taxonomic_sea_all_results.csv",
        sig=f"{IBD}/ibd_taxonomic_sea_results.csv",
    log: f"{OUT}/logs/ibd_taxonomic_sea.log"
    script: "../scripts/ibd_taxonomic_sea.R"

rule ibd_taxonomic_sea_plot:
    input: rules.ibd_taxonomic_sea.output.sig
    output: f"{IBD}/ibd_taxonomic_sea_plot.png"
    params: title="Taxonomic Set Enrichment Analysis (SEA): IBD vs Healthy", ylabel="Taxonomic Group", figsize=(10, 8)
    log: f"{OUT}/logs/ibd_taxonomic_sea_plot.log"
    script: "../scripts/sea_plot.py"

rule ibd_marker_sea:
    output:
        all=f"{IBD}/ibd_marker_clade_sea_all_results.csv",
        sig=f"{IBD}/ibd_marker_clade_sea_results.csv",
    resources: ehub=1  # ExperimentHub cache is not safe for concurrent writers
    log: f"{OUT}/logs/ibd_marker_sea.log"
    script: "../scripts/ibd_marker_sea.R"

rule ibd_marker_sea_plot:
    input: rules.ibd_marker_sea.output.all
    output: f"{IBD}/ibd_marker_sea_plot.png"
    params: title="Marker Clade Enrichment Analysis: IBD vs Healthy", ylabel="Marker Clade", figsize=(12, 10)
    log: f"{OUT}/logs/ibd_marker_sea_plot.log"
    script: "../scripts/sea_plot.py"

rule ibd_coprococcus_extract:
    input: rules.ibd_prepare.output.tse
    output:
        genus=f"{IBD}/coprococcus_genus_abundance.csv",
        species=f"{IBD}/coprococcus_species_abundance.csv",
    log: f"{OUT}/logs/ibd_coprococcus_extract.log"
    script: "../scripts/ibd_coprococcus_extract.R"

rule ibd_coprococcus_stats:
    input:
        genus=rules.ibd_coprococcus_extract.output.genus,
        species=rules.ibd_coprococcus_extract.output.species,
    output:
        stats=f"{IBD}/ibd_coprococcus_species_stats.csv",
        plot=f"{IBD}/ibd_coprococcus_summary_plot.png",
    log: f"{OUT}/logs/ibd_coprococcus_stats.log"
    script: "../scripts/ibd_coprococcus_stats.py"
