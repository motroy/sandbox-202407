"""Alpha diversity (Shannon, Simpson, observed richness) by condition with Mann-Whitney tests."""
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
from ibd_lib import *

setup_log(snakemake)
cfg = snakemake.config
order = list(cfg["condition_labels"].values())      # [control, case]
metrics = cfg["analysis"]["alpha_metrics"]

mat = read_matrix(snakemake.input.matrix)
meta = pd.read_csv(snakemake.input.meta, sep="\t").set_index("sample_id").loc[mat.columns]
p = (mat / mat.sum(axis=0)).to_numpy().T      # samples x features, proportions
with np.errstate(divide="ignore", invalid="ignore"):
    shannon = -np.nansum(np.where(p > 0, p * np.log(p), 0.0), axis=1)
res = pd.DataFrame({"shannon": shannon, "simpson": 1 - (p ** 2).sum(axis=1),
                    "richness": (p > 0).sum(axis=1)}, index=mat.columns)[metrics]
res["condition"] = meta["condition"]
res["dataset"] = meta["dataset"]
res.index.name = "sample_id"
res.to_csv(snakemake.output.table)

sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, len(metrics), figsize=(5.5 * len(metrics), 6), squeeze=False)
for ax, m in zip(axes[0], metrics):
    sns.boxplot(data=res, x="condition", y=m, order=order, hue="condition", hue_order=order, palette="Set2", ax=ax, legend=False)
    sns.stripplot(data=res, x="condition", y=m, order=order, color="black", alpha=0.25, size=2.5, ax=ax)
    Annotator(ax, [tuple(order)], data=res, x="condition", y=m, order=order).configure(
        test="Mann-Whitney", text_format="star", loc="inside", verbose=0).apply_and_annotate()
    ax.set_title(m.capitalize())
    ax.set_xlabel("")
variant = snakemake.wildcards.variant
fig.suptitle(f"Alpha diversity ({'MMUPHin-corrected' if variant == 'mmuphin' else 'uncorrected'})", y=1.02)
plt.tight_layout()
fig.savefig(snakemake.output.plot, dpi=200, bbox_inches="tight")
