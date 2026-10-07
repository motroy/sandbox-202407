"""Raw vs MMUPHin-corrected: variance explained, PCoA, and differential-abundance concordance."""
import sys
sys.path.insert(0, snakemake.scriptdir)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import spearmanr
from ibd_lib import *

setup_log(snakemake)
cfg = snakemake.config
order = list(cfg["condition_labels"].values())
batch_col = cfg["batch_correction"]["batch_column"]
layer = snakemake.wildcards.layer
sns.set_theme(style="whitegrid")

# ---- variance explained (marginal PERMANOVA R2)
rows = []
for stage, f in [("Before", snakemake.input.raw_perm), ("After MMUPHin", snakemake.input.cor_perm)]:
    p = pd.read_csv(f)
    p = p[p.model.str.startswith("batch")] if (p.model.str.startswith("batch")).any() else p
    for _, r in p.iterrows():
        rows.append({"stage": stage, "term": r.term, "R2": r.R2, "p": r.p})
v = pd.DataFrame(rows)
v["term"] = v["term"].replace({batch_col: f"Batch ({batch_col})", "condition": "Condition"})
fig, ax = plt.subplots(figsize=(7, 5))
sns.barplot(data=v, x="term", y="R2", hue="stage", palette=["#999999", "#2a9d8f"], ax=ax)
for c in ax.containers:
    ax.bar_label(c, fmt="%.3f", fontsize=9)
ax.set_ylabel("PERMANOVA R² (Bray-Curtis, marginal)")
ax.set_xlabel("")
ax.set_title(f"Variance explained before / after MMUPHin ({layer})")
ax.legend(title=None)
plt.tight_layout()
fig.savefig(snakemake.output.r2, dpi=200)
plt.close(fig)

# ---- PCoA before/after
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
for j, (stage, f) in enumerate([("Before", snakemake.input.raw_pcoa), ("After MMUPHin", snakemake.input.cor_pcoa)]):
    d = pd.read_csv(f)
    for i, (hue, pal, ho) in enumerate([("batch", "tab10", None), ("condition", "Set2", order)]):
        ax = axes[i, j]
        sns.scatterplot(data=d, x="PCo1", y="PCo2", hue=hue, hue_order=ho, palette=pal, alpha=0.5, s=18, ax=ax)
        ax.set_title(f"{stage} - coloured by {batch_col if hue == 'batch' else 'condition'}")
        ax.set_xlabel(f"PCo1 ({d.var_PCo1.iloc[0]:.1f}%)")
        ax.set_ylabel(f"PCo2 ({d.var_PCo2.iloc[0]:.1f}%)")
        ax.legend(fontsize=8, title=None)
plt.tight_layout()
fig.savefig(snakemake.output.pcoa, dpi=200)
plt.close(fig)

# ---- DA concordance
a = pd.read_csv(snakemake.input.raw_da).set_index("feature")
b = pd.read_csv(snakemake.input.cor_da).set_index("feature")
common = a.index.intersection(b.index)
x, y = a.loc[common, "t_value"], b.loc[common, "t_value"]
ok = x.notna() & y.notna()
rho = spearmanr(x[ok], y[ok])[0] if ok.sum() > 2 else np.nan
thr = cfg["analysis"]["da"]["fdr_threshold"]
sig_raw, sig_cor = int((a.padj < thr).sum()), int((b.padj < thr).sum())
fig, ax = plt.subplots(figsize=(6.5, 6))
ax.scatter(x, y, s=10, alpha=0.5)
lim = np.nanmax(np.abs(np.r_[x, y])) * 1.05
ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.8)
ax.axhline(0, color="grey", lw=0.5); ax.axvline(0, color="grey", lw=0.5)
ax.set_xlabel("t statistic, uncorrected"); ax.set_ylabel("t statistic, MMUPHin-corrected")
ax.set_title(f"DA concordance ({layer}): Spearman ρ = {rho:.2f}\nsignificant: {sig_raw} raw vs {sig_cor} corrected (FDR<{thr})")
plt.tight_layout()
fig.savefig(snakemake.output.da, dpi=200)

summ = v.pivot(index="term", columns="stage", values="R2").reset_index()
summ["layer"] = layer
summ["da_spearman_rho"] = rho
summ["da_significant_raw"], summ["da_significant_corrected"] = sig_raw, sig_cor
summ.to_csv(snakemake.output.summary, index=False)
print(summ.to_string(index=False))
