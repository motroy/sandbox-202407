#!/usr/bin/env python3
"""Map chromosomal blaCTX-M insertion sites onto the E. coli K-12 0-100 min scale.

Input: results/context_hits.tsv (chromosomal rows) + per-sample GBFF (column `gbff` of samples.tsv; LOCUS names must match contig ids).
For each chromosomal hit, up to n_anchors CDS on each side of the inserted region are mapped to K-12
by gene symbol / locus tag / product. The longest contiguous anchor chain per flank (consecutive gap
<= max_chain_gap_min, >= min_chain_genes genes) fixes the position:
  both flanks valid and junction gap <= max_junction_gap_min -> circular midpoint of the nearest anchors
  junction gap larger -> flank with the longer chain (equal chains -> ND)
  one flank valid -> nearest anchor of that flank;  none -> ND
Writes insertion_sites.tsv and site_histogram.tsv (genomes per 0.1 min). Standard library only.
"""
import csv, re, sys, tomllib
from collections import Counter
from pathlib import Path


def parse_cds(gbff):
    """Return {record_name: [(start, end, keys)]} for CDS features; keys = lowercased symbol / locus_tag / product."""
    recs, rec, feat, in_feat = {}, None, None, False

    def close():
        if feat and feat["key"] == "CDS":
            nums = list(map(int, re.findall(r"\d+", feat["loc"])))
            if nums:
                recs[rec].append((min(nums), max(nums), feat["keys"]))

    for line in open(gbff):
        if line.startswith("LOCUS"):
            close(); feat, in_feat = None, False
            rec = line.split()[1]; recs[rec] = []
        elif line.startswith("FEATURES"):
            in_feat = True
        elif line.startswith("ORIGIN"):
            close(); feat, in_feat = None, False
        elif in_feat:
            if line[5:6].strip():                                  # new feature: key at col 6, location at col 22
                close(); feat = {"key": line[5:21].strip(), "loc": line[21:].strip(), "keys": set(), "q": False}
            elif feat:
                t = line[21:].strip()
                if t.startswith("/"):
                    feat["q"] = True
                    m = re.match(r'/(gene|locus_tag|product)="?([^"]*)"?', t)
                    if m: feat["keys"].add(m.group(2).strip().lower())
                elif not feat["q"]:
                    feat["loc"] += t                               # wrapped location (join/complement)
    close()
    return recs


def circ(a, b):                                            # signed circular difference b-a on a 100 scale
    return (b - a + 50) % 100 - 50


def chain(minutes, max_gap, min_genes):
    """Longest run of consecutive anchors (nearest first) with gaps <= max_gap; ties -> nearest."""
    best, run = [], []
    for m in minutes:
        run = run + [m] if run and abs(circ(run[-1], m)) <= max_gap else [m]
        if len(run) > len(best): best = run
    return best if len(best) >= min_genes else []


def main():
    cfg = tomllib.loads(Path(sys.argv[1] if len(sys.argv) > 1 else "config.toml").read_text())
    sc, inp = cfg["sites"], cfg["input"]
    out = Path(inp["outdir"])

    k12 = [c for cs in parse_cds(sc["k12_gbff"]).values() for c in cs]
    k12len = max(e for _, e, _ in k12)
    kmap = {}
    for s, e, keys in k12:
        for k in keys: kmap.setdefault(k, (s + e) / 2 / k12len * 100)

    gbff = {r["sample"]: r.get("gbff") for r in csv.DictReader(open(inp["samples"]), delimiter="\t")}
    rows, cache = [], {}
    for h in csv.DictReader(open(out / "context_hits.tsv"), delimiter="\t"):
        if h["location"] != "chromosome" or not gbff.get(h["sample"]):
            continue
        cds = cache.setdefault(h["sample"], parse_cds(gbff[h["sample"]])).get(h["contig"], [])
        rs, re_ = int(h["region_start"]), int(h["region_end"])
        up = sorted([c for c in cds if c[1] < rs], key=lambda c: -c[1])[:sc["n_anchors"]]
        dn = sorted([c for c in cds if c[0] > re_], key=lambda c: c[0])[:sc["n_anchors"]]

        def flank(cs):
            names = [next((k for k in c[2] if k in kmap), None) for c in cs]
            return chain([kmap[n] for n in names if n], sc["max_chain_gap_min"], sc["min_chain_genes"])
        cu, cd = flank(up), flank(dn)
        pos, how = None, "ND"
        if cu and cd and abs(circ(cu[0], cd[0])) <= sc["max_junction_gap_min"]:
            pos, how = (cu[0] + circ(cu[0], cd[0]) / 2) % 100, "both_flanks"
        elif cu and cd:
            if len(cu) != len(cd):
                pos, how = (cu if len(cu) > len(cd) else cd)[0], "longer_flank"
        elif cu or cd:
            pos, how = (cu or cd)[0], "single_flank"
        rows.append({"sample": h["sample"], "st": h["st"], "gene": h["gene"],
                     "position_min": "ND" if pos is None else f"{round(pos, 1) % 100:.1f}",
                     "method": how, "n_up_chain": len(cu), "n_down_chain": len(cd),
                     "up_anchor_min": f"{cu[0]:.2f}" if cu else "NA",
                     "down_anchor_min": f"{cd[0]:.2f}" if cd else "NA",
                     "mediator": h["mediator"]})

    with open(out / "insertion_sites.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]) if rows else [], delimiter="\t"); w.writeheader(); w.writerows(rows)
    hist = Counter(r["position_min"] for r in rows if r["position_min"] != "ND")
    with open(out / "site_histogram.tsv", "w") as f:
        f.write("position_min\tn_genomes\n")
        for p, n in sorted(hist.items(), key=lambda kv: float(kv[0])): f.write(f"{p}\t{n}\n")
    print(f"{len(rows)} chromosomal insertions, {sum(r['position_min'] == 'ND' for r in rows)} ND", file=sys.stderr)


if __name__ == "__main__":
    main()
