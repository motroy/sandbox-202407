"""Notebook cells 25-26: Coprococcus genus/species abundance IBD vs Healthy
(Mann-Whitney U, BH-FDR across species), figure and descriptive-statistics table."""
import sys as _sys
if snakemake.log:
    _sys.stderr = _sys.stdout = open(snakemake.log[0], "w")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu
from statannotations.Annotator import Annotator
from statsmodels.stats.multitest import multipletests

pc = snakemake.config["ibd"]["coprococcus_pseudocount"]
genus = pd.read_csv(snakemake.input.genus)
species = pd.read_csv(snakemake.input.species)

# Genus
genus_col = [c for c in genus.columns if c not in ("sample_id", "condition")][0]
genus["log_abundance"] = np.log10(genus[genus_col] + pc)

# Species
melted = species.melt(id_vars=["sample_id", "condition"], var_name="Species", value_name="Abundance")
melted["Species_Clean"] = melted["Species"].apply(lambda x: x.split("|")[-1].replace("s__", ""))
melted["log_abundance"] = np.log10(melted["Abundance"] + pc)

species_list = melted["Species_Clean"].unique()
raw_p = []
for s in species_list:
    sub = melted[melted["Species_Clean"] == s]
    _, p = mannwhitneyu(sub.loc[sub.condition == "Healthy", "Abundance"],
                        sub.loc[sub.condition == "IBD", "Abundance"])
    raw_p.append(p)
_, fdr, _, _ = multipletests(raw_p, method="fdr_bh")
fdr_map = dict(zip(species_list, fdr))

sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(1, 2, figsize=(22, 8))

sns.boxplot(data=genus, x="condition", y="log_abundance", hue="condition", palette="Set2", ax=axes[0], legend=False)
sns.stripplot(data=genus, x="condition", y="log_abundance", color="black", alpha=0.05, ax=axes[0])
Annotator(axes[0], [("Healthy", "IBD")], data=genus, x="condition", y="log_abundance").configure(
    test="Mann-Whitney", text_format="star", loc="inside").apply_and_annotate()
axes[0].set_title("Genus Coprococcus (Log10)")
axes[0].set_ylabel("log10(Relative Abundance + 1e-5)")

sns.boxplot(data=melted, x="Species_Clean", y="log_abundance", hue="condition", palette="Set2", ax=axes[1])
axes[1].set_title("Coprococcus Species (Log10)")
axes[1].set_ylabel("log10(Relative Abundance + 1e-5)")
axes[1].legend(title="Condition", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.setp(axes[1].get_xticklabels(), rotation=45, ha="right")
max_val = melted["log_abundance"].max()
axes[1].set_ylim(melted["log_abundance"].min() - 0.5, max_val + 1.5)
for i, s in enumerate(species_list):
    p = fdr_map[s]
    label = "****" if p < 1e-4 else "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "ns"
    axes[1].text(i, max_val + 0.3, label, ha="center", fontweight="bold", fontsize=12)
plt.tight_layout()
fig.savefig(snakemake.output.plot, dpi=300, bbox_inches="tight")

rows = []
for s in species_list:
    for cond in ["Healthy", "IBD"]:
        v = melted[(melted.Species_Clean == s) & (melted.condition == cond)]["Abundance"]
        rows.append({"Species": s, "Condition": cond, "N_Samples": len(v), "Mean": v.mean(),
                     "Median": v.median(), "Std_Dev": v.std(),
                     "IQR": v.quantile(0.75) - v.quantile(0.25),
                     "Stats_Test": "Mann-Whitney U with Benjamini-Hochberg FDR",
                     "padj_fdr_bh": fdr_map[s]})
pd.DataFrame(rows).to_csv(snakemake.output.stats, index=False)
print("--- Species Statistics (FDR BH Corrected) ---")
for s, p in fdr_map.items():
    print(f"{s}: padj = {p:.2e}")
