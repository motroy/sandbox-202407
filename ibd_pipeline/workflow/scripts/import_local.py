"""Standardise one local dataset: MetaPhlAn4 profiles, HUMAnN pathway abundance, sample metadata."""
import sys
sys.path.insert(0, snakemake.scriptdir)
import pandas as pd
from ibd_lib import *

setup_log(snakemake)
spec = snakemake.params.spec
name = spec["name"]
rx = spec.get("sample_id_regex")  # optional regex stripped from file names / column headers

taxa = read_metaphlan(spec["metaphlan"], rx)
print(f"[{name}] MetaPhlAn: {taxa.shape[0]} species x {taxa.shape[1]} samples")

meta = read_metadata(spec["metadata"], spec.get("sample_id_column", "sample_id"))
meta.insert(1, "dataset", name)

missing_meta = sorted(set(taxa.columns) - set(meta["sample_id"]))
if missing_meta:
    print(f"[{name}] WARNING: {len(missing_meta)} MetaPhlAn samples without metadata (dropped): {missing_meta[:5]}...")
common = [s for s in taxa.columns if s in set(meta["sample_id"])]
if not common:
    raise ValueError(f"[{name}] no sample ids shared between MetaPhlAn files and metadata. "
                     f"MetaPhlAn ids look like {list(taxa.columns)[:3]}, metadata ids like {meta['sample_id'].tolist()[:3]}; "
                     "adjust sample_id_regex")
taxa = taxa[common]
meta = meta.set_index("sample_id", drop=False).loc[common].reset_index(drop=True)
write_matrix(taxa, snakemake.output.taxa)
meta.to_csv(snakemake.output.meta, sep="\t", index=False)

if spec.get("humann_pathabundance"):
    path = read_humann(spec["humann_pathabundance"], rx)
    keep = [s for s in path.columns if s in set(common)]
    print(f"[{name}] HUMAnN: {path.shape[0]} unstratified rows x {path.shape[1]} samples ({len(keep)} with MetaPhlAn+metadata)")
    if not keep:
        raise ValueError(f"[{name}] no HUMAnN sample ids match the MetaPhlAn/metadata sample ids; adjust sample_id_regex")
    write_matrix(path[keep], snakemake.output.pathways)
else:
    print(f"[{name}] no humann_pathabundance given; pathways layer will skip this dataset")
    pd.DataFrame({"feature": []}).to_csv(snakemake.output.pathways, sep="\t", index=False)
