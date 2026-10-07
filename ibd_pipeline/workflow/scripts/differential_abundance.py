"""Per-feature differential abundance: log10(abundance + pseudocount) ~ condition [+ batch]; BH-FDR."""
import sys
sys.path.insert(0, snakemake.scriptdir)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from ibd_lib import *

setup_log(snakemake)
cfg = snakemake.config
da = cfg["analysis"]["da"]
control, case = cfg["condition_labels"]["control"], cfg["condition_labels"]["case"]

mat = read_matrix(snakemake.input.matrix)
meta = pd.read_csv(snakemake.input.meta, sep="\t").set_index("sample_id").loc[mat.columns]
y = np.log10(mat.to_numpy() + da["pseudocount"])                  # features x samples
cond = (meta["condition"] == case).astype(float).to_numpy()

cols = [np.ones(len(cond)), cond]
if da["adjust_batch_in_model"]:
    b = pd.get_dummies(meta[cfg["batch_correction"]["batch_column"]], drop_first=True).astype(float)
    cols += [b[c].to_numpy() for c in b.columns]
X = np.column_stack(cols)
dof = X.shape[0] - X.shape[1]
XtXi = np.linalg.pinv(X.T @ X)
beta = y @ X @ XtXi                                               # features x params (XtXi symmetric)
resid = y - beta @ X.T
sigma2 = (resid ** 2).sum(axis=1) / dof
se = np.sqrt(sigma2 * XtXi[1, 1])
with np.errstate(divide="ignore", invalid="ignore"):
    t = beta[:, 1] / se
pval = 2 * stats.t.sf(np.abs(t), dof)

ok = np.isfinite(pval)
padj = np.full(len(pval), np.nan)
from statsmodels.stats.multitest import multipletests
padj[ok] = multipletests(pval[ok], method="fdr_bh")[1]

is_case = (meta["condition"] == case).to_numpy()
out = pd.DataFrame({
    "feature": mat.index, "coef_log10": beta[:, 1], "t_value": t, "pval": pval, "padj": padj,
    f"mean_{control}": mat.to_numpy()[:, ~is_case].mean(axis=1), f"mean_{case}": mat.to_numpy()[:, is_case].mean(axis=1),
    "prevalence": (mat.to_numpy() > 0).mean(axis=1),
}).sort_values("padj")
out.to_csv(snakemake.output.table, index=False)
n_sig = int((out.padj < da["fdr_threshold"]).sum())
print(f"{n_sig} / {len(out)} features with FDR < {da['fdr_threshold']}")


def label(f, layer):
    return species_key(f).replace("s__", "") if layer == "taxa" else (f if len(f) < 60 else f[:57] + "...")


layer, variant = snakemake.wildcards.layer, snakemake.wildcards.variant
top = out.dropna(subset=["t_value"]).reindex(out["t_value"].abs().sort_values(ascending=False).index).head(da["top_n"]).sort_values("t_value")
top["label"] = [label(f, layer) for f in top.feature]
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(10, max(5, 0.32 * len(top) + 2)))
colors = ["#ff9999" if v > 0 else "#66b3ff" for v in top.t_value]
ax.barh(top["label"], top["t_value"], color=colors, edgecolor=["black" if p < da["fdr_threshold"] else "none" for p in top.padj])
ax.axvline(0, color="black", lw=1)
ax.set_xlabel(f"t statistic (← higher in {control} | higher in {case} →)")
ax.set_title(f"Top {len(top)} differentially abundant {layer} ({'MMUPHin-corrected' if variant == 'mmuphin' else 'uncorrected'})\n"
             f"black outline = FDR < {da['fdr_threshold']}  ({n_sig} significant of {len(out)})")
plt.tight_layout()
fig.savefig(snakemake.output.plot, dpi=200)
