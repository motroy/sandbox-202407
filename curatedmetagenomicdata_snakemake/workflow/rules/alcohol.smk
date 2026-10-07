rule alcohol_study:
    output:
        p_alpha=f"{OUT}/alcohol/alcohol_ggplot_alpha.png",
        p_bray=f"{OUT}/alcohol/alcohol_ggplot_bray.png",
        p_umap=f"{OUT}/alcohol/alcohol_ggplot_umap.png",
        meta=f"{OUT}/alcohol/alcohol_meta_alpha.csv",
        bray=f"{OUT}/alcohol/alcohol_bray_coords.csv",
        umap=f"{OUT}/alcohol/alcohol_umap_coords.csv",
    resources: ehub=1  # ExperimentHub cache is not safe for concurrent writers
    log: f"{OUT}/logs/alcohol_study.log"
    script: "../scripts/alcohol_study.R"

rule alcohol_seaborn:
    input:
        meta=rules.alcohol_study.output.meta,
        bray=rules.alcohol_study.output.bray,
        umap=rules.alcohol_study.output.umap,
    output: f"{OUT}/alcohol/alcohol_seaborn_summary.png"
    log: f"{OUT}/logs/alcohol_seaborn.log"
    script: "../scripts/alcohol_plot.py"
