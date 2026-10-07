"""Bray-Curtis PCoA coloured by condition (with ellipses) and by batch; PERMANOVA in the title."""
import sys
sys.path.insert(0, snakemake.scriptdir)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Ellipse
from ibd_lib import *

setup_log(snakemake)
order = list(snakemake.config["condition_labels"].values())
batch_col = snakemake.config["batch_correction"]["batch_column"]
coords = pd.read_csv(snakemake.input.coords)
perm = pd.read_csv(snakemake.input.permanova)


def ellipse(d, ax, color):
    g = d[["PCo1", "PCo2"]]
    if len(g) < 3:
        return
    vals, vecs = np.linalg.eigh(g.cov().to_numpy())
    o = vals.argsort()[::-1]
    vals, vecs = vals[o], vecs[:, o]
    ax.add_patch(Ellipse(xy=g.mean().to_numpy(), width=2 * np.sqrt(vals[0]), height=2 * np.sqrt(vals[1]),
                         angle=np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0])), color=color, alpha=0.2))


def stat(term):
    r = perm[(perm.term == term) & (perm.model != "condition_only")]
    r = r if len(r) else perm[perm.term == term]
    return f"R²={r.R2.iloc[0]:.3f}, p={r.p.iloc[0]:.3g}" if len(r) else "n/a"


sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(16, 7))
pal = dict(zip(order, sns.color_palette("Set2", len(order))))
sns.scatterplot(data=coords, x="PCo1", y="PCo2", hue="condition", hue_order=order, palette=pal, alpha=0.5, s=18, ax=axes[0])
for c in order:
    ellipse(coords[coords.condition == c], axes[0], pal[c])
axes[0].set_title(f"Coloured by condition (PERMANOVA {stat('condition')})")
sns.scatterplot(data=coords, x="PCo1", y="PCo2", hue="batch", palette="tab10", alpha=0.5, s=18, ax=axes[1])
axes[1].set_title(f"Coloured by {batch_col} (PERMANOVA {stat(batch_col)})")
for ax in axes:
    ax.set_xlabel(f"PCo1 ({coords.var_PCo1.iloc[0]:.1f}%)")
    ax.set_ylabel(f"PCo2 ({coords.var_PCo2.iloc[0]:.1f}%)")
v, l = snakemake.wildcards.variant, snakemake.wildcards.layer
fig.suptitle(f"Bray-Curtis PCoA - {l} - {'MMUPHin-corrected' if v == 'mmuphin' else 'uncorrected'}", y=1.0)
plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=200, bbox_inches="tight")
