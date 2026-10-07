"""Figures for the MMUPHin batch correction: PCoA before/after (by study and by condition)
and variance explained (PERMANOVA R2) by study vs condition."""
import sys as _sys
if snakemake.log:
    _sys.stderr = _sys.stdout = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

coords = pd.read_csv(snakemake.input.coords)
perm = pd.read_csv(snakemake.input.permanova)

sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(2, 2, figsize=(15, 12))
for j, stage in enumerate(["Before", "After"]):
    d = coords[coords.stage == stage]
    v1, v2 = d.varexp1.iloc[0], d.varexp2.iloc[0]
    for i, (hue, pal) in enumerate([("study_name", "tab10"), ("condition", "Set2")]):
        ax = axes[i, j]
        sns.scatterplot(data=d, x="PCo1", y="PCo2", hue=hue, palette=pal, alpha=0.5, s=18, ax=ax)
        ax.set_title(f"{stage} MMUPHin - coloured by {'study' if hue == 'study_name' else 'condition'}")
        ax.set_xlabel(f"PCo1 ({v1}%)")
        ax.set_ylabel(f"PCo2 ({v2}%)")
        ax.legend(fontsize=8, title=None)
plt.tight_layout()
fig.savefig(snakemake.output.pcoa, dpi=200)
plt.close(fig)

perm["stage"] = pd.Categorical(perm["stage"], ["Before", "After"])
perm["term"] = perm["term"].map({"study_name": "Study (batch)", "condition": "Condition (Healthy vs IBD)"})
fig, ax = plt.subplots(figsize=(7, 5))
sns.barplot(data=perm, x="term", y="R2", hue="stage", palette=["#999999", "#2a9d8f"], ax=ax)
for cont in ax.containers:
    ax.bar_label(cont, fmt="%.3f", fontsize=9)
ax.set_ylabel("PERMANOVA R$^2$ (Bray-Curtis, marginal)")
ax.set_xlabel("")
ax.set_title("Variance explained before / after MMUPHin")
ax.legend(title=None)
plt.tight_layout()
fig.savefig(snakemake.output.r2, dpi=200)
