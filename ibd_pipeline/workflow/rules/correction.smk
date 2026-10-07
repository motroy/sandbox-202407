# ---- the two analysis "variants": raw (no correction) and mmuphin (batch-corrected)

rule variant_raw:
    input: OUT + "/merged/{layer}/matrix.tsv"
    output: OUT + "/raw/{layer}/matrix.tsv"
    shell: "cp {input} {output}"

rule variant_mmuphin:
    """MMUPHin adjust_batch: batch = batch_correction.batch_column, covariate = condition."""
    input:
        matrix=OUT + "/merged/{layer}/matrix.tsv",
        meta=OUT + "/merged/{layer}/meta.tsv",
    output:
        matrix=OUT + "/mmuphin/{layer}/matrix.tsv",
        info=OUT + "/mmuphin/{layer}/mmuphin_features.csv",
    log: OUT + "/logs/mmuphin_{layer}.log"
    script: "../scripts/mmuphin.R"
