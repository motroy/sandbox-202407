rule demo_asnicar:
    output:
        pca=f"{OUT}/demo/pca.csv",
        alpha=f"{OUT}/demo/alpha_shannon.csv",
        stats=f"{OUT}/demo/demo_stats.txt",
    resources: ehub=1  # ExperimentHub cache is not safe for concurrent writers
    log: f"{OUT}/logs/demo_asnicar.log"
    script: "../scripts/demo_asnicar.R"

rule demo_plot:
    input:
        pca=rules.demo_asnicar.output.pca,
        alpha=rules.demo_asnicar.output.alpha,
    output: f"{OUT}/demo/demo_pca_alpha.png"
    log: f"{OUT}/logs/demo_plot.log"
    script: "../scripts/demo_plot.py"
