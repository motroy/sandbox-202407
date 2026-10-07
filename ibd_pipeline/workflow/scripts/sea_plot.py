"""NES bar plot of significant taxonomic-group enrichment (fgsea)."""
import sys
sys.path.insert(0, snakemake.scriptdir)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from ibd_lib import *

setup_log(snakemake)
cfg = snakemake.config
control, case = cfg["condition_labels"]["control"], cfg["condition_labels"]["case"]
thr = cfg["analysis"]["da"]["fdr_threshold"]
df = pd.read_csv(snakemake.input[0])
df = df[df["padj"] < thr].sort_values("NES") if len(df) else df
fig = plt.figure(figsize=(10, max(4, 0.35 * len(df) + 2)))
v = snakemake.wildcards.variant
title = f"Taxonomic set enrichment: {case} vs {control} ({'MMUPHin-corrected' if v == 'mmuphin' else 'uncorrected'})"
if df.empty:
    plt.text(0.5, 0.5, f"No significant enrichment (FDR < {thr})", ha="center", va="center", fontsize=14)
    plt.axis("off")
else:
    sns.set_theme(style="whitegrid")
    colors = ["#ff9999" if x > 0 else "#66b3ff" for x in df["NES"]]
    sns.barplot(data=df, x="NES", y="pathway", palette=colors, hue="pathway", legend=False)
    plt.xlabel("Normalized Enrichment Score (NES)")
    plt.ylabel("Taxonomic group")
    plt.axvline(0, color="black", lw=1)
    plt.text(df["NES"].min(), -0.8, f"← Enriched in {control}", ha="left", fontweight="bold", color="#004c99")
    plt.text(df["NES"].max(), -0.8, f"Enriched in {case} →", ha="right", fontweight="bold", color="#990000")
plt.title(title, fontsize=13)
plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=200)
