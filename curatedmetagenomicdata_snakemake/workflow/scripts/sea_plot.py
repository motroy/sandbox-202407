"""Notebook cells 21 / 24: horizontal NES bar plot of significant fgsea results."""
import sys as _sys
if snakemake.log:
    _sys.stderr = _sys.stdout = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

df = pd.read_csv(snakemake.input[0])
df = df[df["padj"] < snakemake.config["ibd"]["fdr_threshold"]].sort_values("NES") if "padj" in df and len(df) else df

fig = plt.figure(figsize=tuple(snakemake.params.figsize))
if df.empty:
    plt.text(0.5, 0.5, "No significant enrichment (FDR < %s)" % snakemake.config["ibd"]["fdr_threshold"],
             ha="center", va="center", fontsize=14)
    plt.axis("off")
    plt.title(snakemake.params.title, fontsize=15)
else:
    sns.set_theme(style="whitegrid")
    colors = ["#ff9999" if x > 0 else "#66b3ff" for x in df["NES"]]
    sns.barplot(data=df, x="NES", y="pathway", palette=colors, hue="pathway", legend=False)
    plt.title(snakemake.params.title, fontsize=15)
    plt.xlabel("Normalized Enrichment Score (NES)", fontsize=12)
    plt.ylabel(snakemake.params.ylabel, fontsize=12)
    plt.axvline(0, color="black", linewidth=1)
    plt.text(df["NES"].min(), -0.5, "← Enriched in Healthy", ha="left", va="center", fontweight="bold", color="#004c99")
    plt.text(df["NES"].max(), -0.5, "Enriched in IBD →", ha="right", va="center", fontweight="bold", color="#990000")
plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=300)
