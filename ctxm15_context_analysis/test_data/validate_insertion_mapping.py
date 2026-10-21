#!/usr/bin/env python3
"""Positive control for 03_insertion_sites.py: how well do flank anchors recover a known K-12 position?

Usage: validate_insertion_mapping.py <config.toml> <n_per_genome=100> <sample> [<sample> ...]
For each genome, picks random CDS whose protein has a unique K-12 homolog (its true minute), treats a 3-CDS block centred on
it as an "insertion region", maps the flanks exactly as step 03 does, and compares the placed position with the truth.
Backbone genes only (islands absent from K-12 cannot be scored). Needs blastp/makeblastdb on PATH.
"""
import importlib.util, random, statistics, sys, tempfile, tomllib
from pathlib import Path

spec = importlib.util.spec_from_file_location("sites", Path(__file__).parent.parent / "scripts" / "03_insertion_sites.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
cfg = tomllib.loads(Path(sys.argv[1]).read_text()); sc = cfg["sites"]
sc.setdefault("max_scan", 80); sc.setdefault("min_identity", 80.0); sc.setdefault("min_coverage", 80.0); sc.setdefault("ambiguity_ratio", 1.1)
n = int(sys.argv[2]); samples = sys.argv[3:]
gb = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(cfg["input"]["samples"])}
hdr = gb["sample"]; gcol = hdr.index("gbff")
k12, glen, kname, kid = m.build_k12(sc["k12_gbff"])
random.seed(1)
errs, methods = [], []
for s in samples:
    recs, circ = m.parse_gbff(gb[s][gcol])
    contig = max(recs, key=lambda k: len(recs[k])); cl = sorted(recs[contig], key=lambda c: c["start"])
    with tempfile.TemporaryDirectory() as tmp:
        q = {f"{s}|{c['id']}": c["prot"] for c in cl if c["prot"] and not m.MOBILE.search(c["product"])}
        mapped = m.blast_map(q, k12, kid, sc, tmp)
    def minute(c):
        return mapped.get(f"{s}|{c['id']}")
    idx = [i for i in range(2, len(cl) - 2) if minute(cl[i]) is not None]
    for i in random.sample(idx, min(n, len(idx))):
        truth = minute(cl[i]); rs, re_ = cl[i - 1]["start"], cl[i + 1]["end"]
        up, dn = m.flank_candidates(cl, rs, re_, contig in circ, sc["max_scan"])
        def anchors(side):
            ms = []
            for c in side:
                if m.MOBILE.search(c["product"]): continue
                v = minute(c)
                if v is not None: ms.append(v)
                if len(ms) >= sc["n_anchors"]: break
            return ms
        pos, how, _, _ = m.place(anchors(up), anchors(dn), sc)
        methods.append(how)
        if pos is not None: errs.append(abs(m.circ(truth, pos)))
tot = len(methods)
print(f"{tot} pseudo-insertions in {len(samples)} genomes; placed {len(errs)} ({100 * len(errs) / tot:.1f}%), ND {tot - len(errs)}")
if errs:
    errs.sort()
    print(f"|error| median {statistics.median(errs):.2f} min, 90th pct {errs[int(0.9 * len(errs)) - 1]:.2f}, max {errs[-1]:.2f}; "
          f"within 0.5 min {100 * sum(e <= 0.5 for e in errs) / len(errs):.1f}%, within 1 min {100 * sum(e <= 1 for e in errs) / len(errs):.1f}%")
from collections import Counter
print(dict(Counter(methods)))
