"""Define case/control, filter samples, join datasets and write one matrix + metadata per layer."""
import sys
sys.path.insert(0, snakemake.scriptdir)
import numpy as np
import pandas as pd
from ibd_lib import *

setup_log(snakemake)
cfg = snakemake.config
labels = cfg["condition_labels"]
join, renorm = cfg["merge"]["join"], cfg["merge"]["renormalize"]
drop_prefixes = tuple(cfg["merge"].get("humann_drop") or [])
cur = cfg["datasets"]["curated"]
datasets, curated, local = snakemake.params.datasets, set(snakemake.params.curated), snakemake.params.local

taxa_in = dict(zip(datasets, snakemake.input.taxa))
path_in = dict(zip(datasets, snakemake.input.pathways))
meta_in = dict(zip(datasets, snakemake.input.meta))

# ---- per-dataset sample table with the unified 'condition'
metas = {}
for ds in datasets:
    m = pd.read_csv(meta_in[ds], sep="\t", dtype=str, keep_default_na=True)
    spec = cur if ds in curated else local[ds]
    col = spec["condition_column"]
    if col not in m.columns:
        raise KeyError(f"[{ds}] condition_column '{col}' not in metadata columns")
    if ds in curated and cur.get("body_site") and "body_site" in m.columns:
        m = m[m["body_site"] == cur["body_site"]]
    m = m.copy()
    m["condition"] = map_condition(m[col], spec, labels)
    n_drop = int(m["condition"].isna().sum())
    print(f"[{ds}] {len(m)} samples; {n_drop} dropped (not matching case/control definition; values: "
          f"{dict(m.loc[m['condition'].isna(), col].value_counts(dropna=False).head(5))})")
    metas[ds] = m[m["condition"].notna()].set_index("sample_id", drop=False)

all_ids = pd.concat([m["sample_id"] for m in metas.values()])
if all_ids.duplicated().any():
    raise ValueError(f"Sample ids are not unique across datasets, e.g. {all_ids[all_ids.duplicated()].tolist()[:5]}")

batch_col = cfg["batch_correction"]["batch_column"]


def build_layer(layer):
    tables = {}
    for ds in datasets:
        path = (taxa_in if layer == "taxa" else path_in)[ds]
        df = pd.read_csv(path, sep="\t", index_col=0)
        if df.shape[1] == 0 or df.shape[0] == 0:
            print(f"[{layer}] {ds}: no data, skipped")
            continue
        df = df.loc[:, [s for s in df.columns if s in metas[ds].index]]
        if layer == "pathways" and drop_prefixes:
            df = df.loc[~df.index.str.startswith(drop_prefixes)]
        if df.shape[1] == 0:
            print(f"[{layer}] {ds}: no samples left after filtering, skipped")
            continue
        tables[ds] = df
    if not tables:
        raise ValueError(f"Layer '{layer}': none of the selected datasets provides data")

    # taxa are matched on the species name (robust to different MetaPhlAn database lineages)
    lineage = {}
    keyed = {}
    for ds, df in tables.items():
        if layer == "taxa":
            keys = [species_key(i) for i in df.index]
            for k, i in zip(keys, df.index):
                lineage.setdefault(k, i)
            df = df.copy(); df.index = keys
            df = df.groupby(level=0).sum()
        keyed[ds] = df
    if join == "inner":
        feats = sorted(set.intersection(*[set(d.index) for d in keyed.values()]))
    else:
        feats = sorted(set.union(*[set(d.index) for d in keyed.values()]))
    if not feats:
        raise ValueError(f"Layer '{layer}': no shared features between datasets (try merge.join: outer)")
    mat = pd.concat([d.reindex(feats).fillna(0.0) for d in keyed.values()], axis=1)
    if layer == "taxa":
        mat.index = [lineage[k] for k in mat.index]
    mat = mat.loc[:, mat.sum(axis=0) > 0]            # drop empty samples
    if renorm:
        mat = mat / mat.sum(axis=0) * 100.0
    meta = pd.concat([metas[ds] for ds in tables], axis=0).loc[mat.columns]
    meta.index.name = None
    if batch_col not in meta.columns or meta[batch_col].isna().any():
        raise KeyError(f"batch_correction.batch_column '{batch_col}' is missing/NA for some samples in layer '{layer}'")
    print(f"[{layer}] datasets: {list(tables)}; {mat.shape[0]} features x {mat.shape[1]} samples ({join} join)")
    return mat, meta


counts = []
for layer, out_m, out_meta in [("taxa", snakemake.output.taxa, snakemake.output.taxa_meta),
                               ("pathways", snakemake.output.pathways, snakemake.output.pathways_meta)]:
    try:
        mat, meta = build_layer(layer)
    except ValueError as e:
        if layer == "pathways" and "none of the selected datasets" in str(e) and "pathways" not in snakemake.config["layers"]:
            print(f"[pathways] {e} (layer not requested; writing empty placeholder)")
            pd.DataFrame({"feature": []}).to_csv(out_m, sep="\t", index=False)
            pd.DataFrame({"sample_id": []}).to_csv(out_meta, sep="\t", index=False)
            continue
        raise
    write_matrix(mat, out_m)
    meta.to_csv(out_meta, sep="\t", index=False)
    c = meta.groupby(["dataset", "condition"]).size().rename("n_samples").reset_index()
    c.insert(0, "layer", layer)
    counts.append(c)
pd.concat(counts).to_csv(snakemake.output.counts, index=False)
print(pd.concat(counts).to_string(index=False))
