#!/usr/bin/env python3
"""Stratified test subset of the Zenodo 20934804 table (2,507 ESBL E. coli genomes; run -I for safety).

Usage: make_subset.py summary_df_2507genomes.txt subset.tsv [seed]
Strata (n genomes): chromosomal blaCTX-M ST131 (5), ST38 (4), other STs (5), chromosomal with a multi-ARG
barcode (4), plasmid-only ST131 (4), plasmid-only other (4), both chromosomal and plasmidic (4). Seeded; no duplicates.
"""
import csv, random, sys

src, out = sys.argv[1], sys.argv[2]
random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
rows = list(csv.DictReader(open(src), delimiter="\t"))
none = lambda v: v in ("none", "", "NA")
chrom = lambda r: not none(r["Chromosomal CTX-M variant"])
plas = lambda r: not none(r["Plasmidic CTX-M variant"])
multi = lambda r: r["ARGs barcode"].count("]_") >= 1
strata = [
    ("chrom_ST131", lambda r: chrom(r) and r["Sequence Type"] == "131", 5),
    ("chrom_ST38", lambda r: chrom(r) and r["Sequence Type"] == "38", 4),
    ("chrom_other", lambda r: chrom(r) and r["Sequence Type"] not in ("131", "38", "ND"), 5),
    ("chrom_multiARG", lambda r: chrom(r) and multi(r), 4),
    ("plasmid_ST131", lambda r: plas(r) and not chrom(r) and r["Sequence Type"] == "131", 4),
    ("plasmid_other", lambda r: plas(r) and not chrom(r) and r["Sequence Type"] not in ("131", "ND"), 4),
    ("both", lambda r: chrom(r) and plas(r), 4),
]
picked, seen = [], set()
for name, f, n in strata:
    pool = [r for r in rows if f(r) and r["Strain"] not in seen]
    for r in random.sample(pool, n):
        seen.add(r["Strain"]); picked.append({"stratum": name, **r})
with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, list(picked[0]), delimiter="\t"); w.writeheader(); w.writerows(picked)
print(f"{len(picked)} genomes -> {out}", file=sys.stderr)
