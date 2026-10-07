"""Notebook cell 19: skbio PERMANOVA (Healthy vs IBD, Bray-Curtis) + within/between distance diagnostics."""
import sys as _sys
if snakemake.log:
    _sys.stderr = _sys.stdout = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from skbio.stats.distance import DistanceMatrix, permanova
from statannotations.Annotator import Annotator

import sys
sys.path.insert(0, snakemake.scriptdir)
from ibd_alpha_beta_plots_helpers import draw_ellipse  # noqa: E402

cfg = snakemake.config["ibd"]
np.random.seed(snakemake.config["seed"])

plot_df = pd.read_csv(snakemake.input.plot_df, index_col=0)
pcoa = pd.read_csv(snakemake.input.pcoa, index_col=0)
dist = pd.read_csv(snakemake.input.dist, index_col=0)

# Restrict distance matrix to the samples used in the plots (same order)
ids = [i for i in plot_df.index if i in dist.index]
plot_df = plot_df.loc[ids]
mat = dist.loc[ids, ids].to_numpy(dtype=float)
pcoa_plot = pcoa.loc[ids].copy()
pcoa_plot["condition"] = plot_df["condition"]

dm = DistanceMatrix(mat, ids=ids)
res = permanova(dm, plot_df[["condition"]], column="condition", permutations=cfg["permutations"])
p_val = res["p-value"]
with open(snakemake.output.txt, "w") as fh:
    fh.write("PERMANOVA (scikit-bio), Bray-Curtis ~ condition (Healthy vs IBD)\n")
    fh.write(res.to_string() + "\n")

pos = {s: i for i, s in enumerate(ids)}
h = [pos[s] for s in plot_df.index[plot_df["condition"] == "Healthy"]]
i_ = [pos[s] for s in plot_df.index[plot_df["condition"] == "IBD"]]
within_h = mat[np.ix_(h, h)][np.triu_indices(len(h), k=1)]
within_i = mat[np.ix_(i_, i_)][np.triu_indices(len(i_), k=1)]
between = mat[np.ix_(h, i_)].flatten()
dist_df = pd.DataFrame({
    "Distance": np.concatenate([within_h, within_i, between]),
    "Comparison": ["Within Healthy"] * len(within_h) + ["Within IBD"] * len(within_i) + ["Between Groups"] * len(between),
})

sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(18, 7))
sns.scatterplot(data=pcoa_plot, x=pcoa_plot.columns[0], y=pcoa_plot.columns[1], hue="condition",
                palette="Set2", alpha=0.4, ax=axes[0])
for cond, color in {"Healthy": "orange", "IBD": "green"}.items():
    if cond in pcoa_plot["condition"].values:
        draw_ellipse(pcoa_plot, cond, axes[0], color)
axes[0].set_title(f"Bray-Curtis PCoA (PERMANOVA p = {p_val:.4f})")

sns.boxplot(data=dist_df, x="Comparison", y="Distance", hue="Comparison", palette="Set2", ax=axes[1], legend=False)
pairs = [("Within Healthy", "Between Groups"), ("Within IBD", "Between Groups"), ("Within Healthy", "Within IBD")]
Annotator(axes[1], pairs, data=dist_df, x="Comparison", y="Distance").configure(
    test="Mann-Whitney", text_format="star", loc="inside").apply_and_annotate()
axes[1].set_title("Pairwise Distance Distribution")
plt.tight_layout()
fig.savefig(snakemake.output.fig, dpi=300)
