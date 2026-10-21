#!/usr/bin/env python3
"""Compare pipeline output with the paper's per-genome calls (Zenodo table columns).

Usage: compare_to_paper.py subset.tsv results_dir   -> prints a table and writes comparison.tsv
Compared: n chromosomal / plasmidic blaCTX-M copies, the chromosomal variant set, the chromosomal ARG barcode
(AMRFinderPlus names differ from the paper's ResFinder names, so barcodes are shown, not auto-matched),
ST, and ST131 clade.
"""
import csv, re, sys
from collections import defaultdict
from pathlib import Path

sub = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
res = Path(sys.argv[2])
hits = defaultdict(list)
for h in csv.DictReader(open(res / "context_hits.tsv"), delimiter="\t"):
    hits[h["sample"]].append(h)
typ = {r["sample"]: r for r in csv.DictReader(open(res / "typing.tsv"), delimiter="\t")} if (res / "typing.tsv").exists() else {}

def paper(v):
    d = {}
    for m in re.finditer(r"(blaCTX-M-[\w]+) \(n=(\d+)\)", v): d[m.group(1)] = int(m.group(2))
    return d

out, agree = [], 0
for r in sub:
    s = r["Strain"]; hs = hits.get(s, [])
    pc, pp = paper(r["Chromosomal CTX-M variant"]), paper(r["Plasmidic CTX-M variant"])
    oc = [h for h in hs if h["location"] == "chromosome"]; op = [h for h in hs if h["location"] == "plasmid"]
    ok = (sum(pc.values()) == len(oc)) and (sum(pp.values()) == len(op))
    agree += ok
    t = typ.get(s, {})
    out.append({"sample": s, "stratum": r["stratum"],
                "paper_chrom": sum(pc.values()), "ours_chrom": len(oc), "paper_plasmid": sum(pp.values()),
                "ours_plasmid": len(op), "copies_agree": int(ok),
                "paper_chrom_variants": ";".join(sorted(pc)), "ours_chrom_variants": ";".join(sorted({h["gene"] for h in oc})),
                "paper_st": "ND" if r["Sequence Type"] == "ND" else "ST" + r["Sequence Type"], "ours_st": t.get("st", "NA"),
                "paper_st131_clade": r["Sequence Type 131 clade"], "ours_st131_subclade": t.get("st131_subclade", "NA"),
                "paper_barcode": r["ARGs barcode"], "ours_barcodes": " | ".join(sorted({h["barcode"] for h in oc})) or "none",
                "ours_mediators": ",".join(sorted({h["mediator"] for h in oc})) or "none"})
# Barcode check: genes co-located with chromosomal blaCTX-M (paper: "blaCTX-M-15: {[a]_[b]}" strings; ours: ARGs in the window)
def partners(txt):
    return {g for g in re.findall(r"\[([^\]]+)\]", txt) if not g.startswith("blaCTX-M")}
bc_rows = []
for o, r in zip(out, sub):
    if r["ARGs barcode"] in ("none", "Not evaluated"):
        continue
    ps = partners(r["ARGs barcode"])
    ours_p = {g for h in hits.get(r["Strain"], []) if h["location"] == "chromosome"
              for g in h["co_resistance_genes"].split(";") if g not in ("none",) and not g.startswith("blaCTX-M")}
    j = len(ps & ours_p) / len(ps | ours_p) if (ps | ours_p) else 1.0
    bc_rows.append((r["Strain"], ps, ours_p, j))
    o["barcode_partners_paper"] = ";".join(sorted(ps)); o["barcode_partners_ours"] = ";".join(sorted(ours_p)); o["barcode_jaccard"] = f"{j:.2f}"
if bc_rows:
    exact = sum(j == 1.0 for *_, j in bc_rows)
    print(f"ARG barcode partners (genomes with a paper barcode): exact {exact}/{len(bc_rows)}, mean Jaccard {sum(j for *_, j in bc_rows) / len(bc_rows):.2f}")
    for sid, ps, op, j in bc_rows:
        if j < 1.0: print(f"  {sid}: paper-only {sorted(ps - op)} ours-only {sorted(op - ps)}")
for o in out: o.setdefault("barcode_jaccard", "NA")
with open(res / "comparison.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, list(out[0]), delimiter="\t"); w.writeheader(); w.writerows(out)
print(f"copy numbers (chromosomal and plasmid) agree with the paper in {agree}/{len(sub)} genomes")
for o in out:
    if not o["copies_agree"]:
        print(f"  {o['sample']} ({o['stratum']}): chrom paper={o['paper_chrom']} ours={o['ours_chrom']}; plasmid paper={o['paper_plasmid']} ours={o['ours_plasmid']}")
print(f"ST agree: {sum(o['paper_st'] == o['ours_st'] for o in out)}/{len(out)}; "
      f"ST131 clade agree: {sum(o['paper_st131_clade'] == o['ours_st131_subclade'] for o in out if o['paper_st'] == 'ST131')}/{sum(o['paper_st'] == 'ST131' for o in out)}")
