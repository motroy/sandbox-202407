"""Notebook cell 11: seaborn re-plot of the alcohol study (alpha, Bray-Curtis PCoA, UMAP)."""
import sys as _sys
if snakemake.log:
    _sys.stderr = _sys.stdout = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

meta = pd.read_csv(snakemake.input.meta)
bray = pd.read_csv(snakemake.input.bray).merge(meta[["sample_id", "alcohol"]], on="sample_id")
umap = pd.read_csv(snakemake.input.umap).merge(meta[["sample_id", "alcohol"]], on="sample_id")

sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 3, figsize=(20, 6))

sns.boxplot(data=meta, x="alcohol", y="shannon", hue="alcohol", palette="Set2", ax=axes[0], legend=False)
sns.stripplot(data=meta, x="alcohol", y="shannon", color="black", alpha=0.3, ax=axes[0])
axes[0].set_title("Alpha Diversity (Shannon Index)")

for ax, df, title, lab in [(axes[1], bray, "Bray-Curtis PCoA (Python/Seaborn)", "PCo"),
                           (axes[2], umap, "UMAP (Python/Seaborn)", "UMAP")]:
    xcol, ycol = df.columns[0], df.columns[1]
    sns.scatterplot(data=df, x=xcol, y=ycol, hue="alcohol", style="alcohol", s=80, ax=ax)
    ax.set_title(title)
    ax.set_xlabel(f"{lab} 1")
    ax.set_ylabel(f"{lab} 2")

plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=150)
