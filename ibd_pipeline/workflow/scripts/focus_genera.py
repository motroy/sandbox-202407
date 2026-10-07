"""Focused plots for selected genera: genus total and per-species abundance by condition (Mann-Whitney, BH-FDR)."""
import sys
sys.path.insert(0, snakemake.scriptdir)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests
from ibd_lib import *

setup_log(snakemake)
cfg = snakemake.config
order = list(cfg["condition_labels"].values())
pc = cfg["analysis"]["focus_pseudocount"]
genera = cfg["analysis"]["focus_genera"]

mat = read_matrix(snakemake.input.matrix)
meta = pd.read_csv(snakemake.input.meta, sep="\t").set_index("sample_id").loc[mat.columns]
genus_of = pd.Series({f: parse_lineage(f).get("Genus") for f in mat.index})

sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(max(len(genera), 1), 2, figsize=(18, 7 * max(len(genera), 1)), squeeze=False,
                         gridspec_kw={"width_ratios": [1, 2]})
rows = []
for ax_row, g in zip(axes, genera):
    sp = genus_of.index[genus_of == g]
    if len(sp) == 0:
        for ax in ax_row:
            ax.text(0.5, 0.5, f"Genus {g} not found", ha="center"); ax.axis("off")
        continue
    gen = mat.loc[sp].sum(axis=0)
    gdf = pd.DataFrame({"condition": meta["condition"], "log_abundance": np.log10(gen + pc)})
    sns.boxplot(data=gdf, x="condition", y="log_abundance", order=order, hue="condition", hue_order=order, palette="Set2", ax=ax_row[0], legend=False)
    sns.stripplot(data=gdf, x="condition", y="log_abundance", order=order, color="black", alpha=0.05, ax=ax_row[0])
    Annotator(ax_row[0], [tuple(order)], data=gdf, x="condition", y="log_abundance", order=order).configure(
        test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()
    ax_row[0].set_title(f"Genus {g} (log10)")
    ax_row[0].set_ylabel(f"log10(relative abundance + {pc:g})")

    long = mat.loc[sp].T.assign(condition=meta["condition"].values).melt(id_vars="condition", var_name="feature", value_name="abundance")
    long["Species"] = long["feature"].map(lambda f: species_key(f).replace("s__", ""))
    long["log_abundance"] = np.log10(long["abundance"] + pc)
    pvals, names = [], []
    for s, d in long.groupby("Species"):
        a, b = d.loc[d.condition == order[0], "abundance"], d.loc[d.condition == order[1], "abundance"]
        p = mannwhitneyu(a, b).pvalue if a.nunique() + b.nunique() > 1 else np.nan
        names.append(s); pvals.append(p)
        rows.append({"genus": g, "species": s, **{f"n_{c}": int((d.condition == c).sum()) for c in order},
                     **{f"median_{c}": d.loc[d.condition == c, "abundance"].median() for c in order},
                     **{f"mean_{c}": d.loc[d.condition == c, "abundance"].mean() for c in order}, "pval": p})
    pv = np.array(pvals, float)
    padj = np.full(len(pv), np.nan)
    ok = np.isfinite(pv)
    if ok.any():
        padj[ok] = multipletests(pv[ok], method="fdr_bh")[1]
    for r, pa in zip(rows[-len(names):], padj):
        r["padj_fdr_bh"] = pa
    sns.boxplot(data=long, x="Species", y="log_abundance", hue="condition", hue_order=order, palette="Set2", ax=ax_row[1])
    ax_row[1].set_title(f"{g} species (log10; stars = BH-FDR)")
    ax_row[1].set_ylabel(f"log10(relative abundance + {pc:g})")
    plt.setp(ax_row[1].get_xticklabels(), rotation=45, ha="right")
    ymax = long["log_abundance"].max()
    ax_row[1].set_ylim(long["log_abundance"].min() - 0.5, ymax + 1.5)
    for i, (s, pa) in enumerate(zip(names, padj)):
        lab = "ns" if not np.isfinite(pa) else "****" if pa < 1e-4 else "***" if pa < 1e-3 else "**" if pa < 1e-2 else "*" if pa < 0.05 else "ns"
        ax_row[1].text(i, ymax + 0.3, lab, ha="center", fontweight="bold")
plt.tight_layout()
fig.savefig(snakemake.output.plot, dpi=200, bbox_inches="tight")
pd.DataFrame(rows).to_csv(snakemake.output.stats, index=False)
