"""Notebook cell 7: PCA + Shannon-by-body-site figure for the AsnicarF_2017 demo."""
import sys as _sys
if snakemake.log:
    _sys.stderr = _sys.stdout = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

pca = pd.read_csv(snakemake.input.pca)
alpha = pd.read_csv(snakemake.input.alpha)

sns.set_theme(style="whitegrid")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

sns.scatterplot(data=pca, x="PC1", y="PC2", hue="body_site", style="body_site", s=100, ax=ax1)
ax1.set_title("PCA of Microbial Relative Abundance")

sns.boxplot(data=alpha, x="body_site", y="shannon", hue="body_site", palette="Set2", ax=ax2, legend=False)
sns.stripplot(data=alpha, x="body_site", y="shannon", color="black", alpha=0.5, ax=ax2)
ax2.set_title("Shannon Diversity by Body Site")
ax2.set_xlabel("Body Site")
ax2.set_ylabel("Shannon Index")

plt.tight_layout()
fig.savefig(snakemake.output[0], dpi=150)
