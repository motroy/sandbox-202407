# ---- standardise every dataset to  taxa.tsv / pathways.tsv / meta.tsv  (features x samples, percent)

rule fetch_curated:
    """Download one curatedMetagenomicData study (relative_abundance + pathway_abundance)."""
    output:
        taxa=DATA + "/curated/{study}/taxa.tsv",
        pathways=DATA + "/curated/{study}/pathways.tsv",
        meta=DATA + "/curated/{study}/meta.tsv",
    wildcard_constraints: study="|".join(CURATED) if CURATED else "NONE"
    resources: ehub=1   # ExperimentHub cache is not safe for concurrent writers
    log: OUT + "/logs/fetch_curated_{study}.log"
    script: "../scripts/fetch_curated.R"

rule import_local:
    """Parse local MetaPhlAn4 profiles and HUMAnN pathway abundance + sample metadata."""
    input: lambda wc: local_inputs(wc.name)
    output:
        taxa=DATA + "/local/{name}/taxa.tsv",
        pathways=DATA + "/local/{name}/pathways.tsv",
        meta=DATA + "/local/{name}/meta.tsv",
    wildcard_constraints: name="|".join(LOCAL) if LOCAL else "NONE"
    params: spec=lambda wc: LOCAL[wc.name]
    log: OUT + "/logs/import_local_{name}.log"
    script: "../scripts/import_local.py"

rule merge_datasets:
    """Apply case/control definitions, join datasets, build per-layer matrix + metadata."""
    input:
        taxa=[f"{dataset_dir(d)}/taxa.tsv" for d in DATASETS],
        pathways=[f"{dataset_dir(d)}/pathways.tsv" for d in DATASETS],
        meta=[f"{dataset_dir(d)}/meta.tsv" for d in DATASETS],
    output:
        counts=OUT + "/merged/sample_counts.csv",
        taxa=OUT + "/merged/taxa/matrix.tsv",
        taxa_meta=OUT + "/merged/taxa/meta.tsv",
        pathways=OUT + "/merged/pathways/matrix.tsv",
        pathways_meta=OUT + "/merged/pathways/meta.tsv",
    params:
        datasets=DATASETS,
        curated=CURATED,
        local=LOCAL,
    log: OUT + "/logs/merge_datasets.log"
    script: "../scripts/merge_datasets.py"
