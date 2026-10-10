"""Compare Amira (reads) and AMRFinderPlus (assembly) AMR gene calls at gene-family level.

usage: compare.py AMIRA_RESULTS.tsv AMRFINDER.tsv OUTDIR
"""
import re
import sys
from pathlib import Path

import pandas as pd

amira_tsv, amrf_tsv, outdir = sys.argv[1:4]
out = Path(outdir)
out.mkdir(parents=True, exist_ok=True)


def family(name: str) -> str:
    """Collapse allele-level names to a gene family so the tools' naming is comparable."""
    n = name.strip()
    n = re.sub(r"^bla(CTX-M|TEM|OXA|SHV|KPC|CMY|NDM|EC|VIM|IMP)-?.*$", r"bla\1", n)  # blaTEM-1 -> blaTEM
    n = re.sub(r"-cr\d*$", "", n)  # aac(6')-Ib-cr5 -> aac(6')-Ib
    n = re.sub(r"^(aac\(6'\)-Ib)\d+$", r"\1", n)  # aac(6')-Ib4 -> aac(6')-Ib
    n = re.sub(r"^(dfrA)(\d+)$", r"\1\2", n)
    return n


amira = pd.read_csv(amira_tsv, sep="\t")
amrf = pd.read_csv(amrf_tsv, sep="\t")

# AMRFinderPlus: acquired AMR genes only (point mutations and non-AMR classes are out of Amira's scope)
amrf_amr = amrf[(amrf["Type"] == "AMR") & (amrf["Subtype"] == "AMR")].copy()
amrf_amr["family"] = amrf_amr["Element symbol"].map(family)

# Amira: one row per allele; a row can carry several equally-good reference names
rows = []
for _, r in amira.iterrows():
    for nm in str(r["Determinant name"]).split("/"):
        rows.append({"family": family(nm), "allele": r["Amira allele"], "depth": r["Relative mean read depth"]})
amira_long = pd.DataFrame(rows).drop_duplicates(["family", "allele"])

fam_amira = amira_long.groupby("family").agg(amira_alleles=("allele", "nunique"), amira_sum_depth=("depth", "sum"))
fam_amrf = amrf_amr.groupby("family").agg(
    amrfinder_copies=("Contig id", "size"),
    amrfinder_symbols=("Element symbol", lambda s: ";".join(sorted(set(s)))),
    amrfinder_scope=("Scope", lambda s: ";".join(sorted(set(s)))),
)
cmp_ = fam_amira.join(fam_amrf, how="outer")
cmp_["amira_alleles"] = cmp_["amira_alleles"].fillna(0).astype(int)
cmp_["amrfinder_copies"] = cmp_["amrfinder_copies"].fillna(0).astype(int)


def status(r):
    if r.amira_alleles and r.amrfinder_copies:
        return "both"
    return "amira_only" if r.amira_alleles else "amrfinderplus_only"


cmp_["status"] = cmp_.apply(status, axis=1)
cmp_ = cmp_.sort_values(["status", "family"])
cmp_.round(3).to_csv(out / "gene_family_comparison.tsv", sep="\t")

n = cmp_["status"].value_counts().to_dict()
with open(out / "summary.txt", "w") as fh:
    fh.write(f"Amira allele rows: {len(amira)}; distinct gene families: {len(fam_amira)}\n")
    fh.write(f"AMRFinderPlus AMR (acquired, Subtype=AMR) hits: {len(amrf_amr)}; families: {len(fam_amrf)}\n")
    fh.write(f"Shared: {n.get('both', 0)}  Amira-only: {n.get('amira_only', 0)}  AMRFinderPlus-only: {n.get('amrfinderplus_only', 0)}\n")
print(open(out / "summary.txt").read())
print(cmp_.round(2).to_string())
