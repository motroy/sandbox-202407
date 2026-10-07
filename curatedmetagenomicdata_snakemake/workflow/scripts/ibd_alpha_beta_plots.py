"""Notebook cells 17-18: IBD alpha diversity (Shannon, Faith PD) and Bray-Curtis PCoA with ellipses."""
import sys as _sys
if snakemake.log:
    _sys.stderr = _sys.stdout = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from statannotations.Annotator import Annotator

import sys
sys.path.insert(0, snakemake.scriptdir)
from ibd_alpha_beta_plots_helpers import draw_ellipse  # noqa: E402

meta = pd.read_csv(snakemake.input.meta, index_col=0)
pcoa = pd.read_csv(snakemake.input.pcoa, index_col=0)

diversity_cols = [c for c in ["shannon_index", "faiths_pd"] if c in meta.columns]
plot_df = meta.dropna(subset=diversity_cols + ["condition"]).copy()
plot_df.index.name = "sample_id"
plot_df.to_csv(snakemake.output.plot_df)
faith_col = "faiths_pd" if "faiths_pd" in plot_df.columns else None

PAIR = [("Healthy", "IBD")]


def alpha_panel(ax, col, title):
    sns.boxplot(data=plot_df, x="condition", y=col, hue="condition", palette="Set2", ax=ax, legend=False)
    sns.stripplot(data=plot_df, x="condition", y=col, color="black", alpha=0.3, ax=ax)
    Annotator(ax, PAIR, data=plot_df, x="condition", y=col).configure(
        test="Mann-Whitney", text_format="star", loc="inside").apply_and_annotate()
    ax.set_title(title)


sns.set_theme(style="whitegrid")

# Overview: Shannon | Faith PD | PCoA
fig, axes = plt.subplots(1, 3, figsize=(22, 6))
alpha_panel(axes[0], "shannon_index", "Shannon Diversity")
if faith_col:
    alpha_panel(axes[1], faith_col, "Faith's Phylogenetic Diversity")
else:
    axes[1].text(0.5, 0.5, "Faith's PD column not found", ha="center")

pcoa_plot = pcoa.loc[plot_df.index].copy()
pcoa_plot["condition"] = plot_df["condition"]
sns.scatterplot(data=pcoa_plot, x=pcoa_plot.columns[0], y=pcoa_plot.columns[1], hue="condition",
                palette="Set2", alpha=0.4, ax=axes[2])
for cond, color in {"Healthy": "orange", "IBD": "green"}.items():
    if cond in pcoa_plot["condition"].values:
        draw_ellipse(pcoa_plot, cond, axes[2], color)
axes[2].set_title("Bray-Curtis PCoA")
plt.tight_layout()
fig.savefig(snakemake.output.overview, dpi=150, bbox_inches="tight")
plt.close(fig)

# Alpha-only figure (cell 18)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(15, 6))
alpha_panel(a1, "shannon_index", "Shannon Diversity")
if faith_col:
    alpha_panel(a2, faith_col, "Faith's Phylogenetic Diversity")
plt.tight_layout()
fig.savefig(snakemake.output.alpha, dpi=300, bbox_inches="tight")
