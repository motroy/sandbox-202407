#!/usr/bin/env python3
"""Run ISEScan on windows around each blaCTX-M hit (whole-genome ISEScan takes 15+ min per 5 Mb genome here).

Usage: run_isescan_windows.py <results/context_hits.tsv> <genomes_dir> <outdir> [flank_bp=30000]
Per genome: one FASTA with a window per hit (header "<contig>@<offset>", original coordinate = local + offset),
ISEScan on it -> <outdir>/<sample>/<sample>.windows.fa.tsv. Needs isescan.py on PATH.
"""
import csv, subprocess, sys
from collections import defaultdict
from pathlib import Path

hits_file, gdir, out = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
flank = int(sys.argv[4]) if len(sys.argv) > 4 else 30000
out = out.resolve(); out.mkdir(parents=True, exist_ok=True); gdir = gdir.resolve()


def read_fa(p):
    seqs, n = {}, None
    for line in open(p):
        if line.startswith(">"): n = line[1:].split()[0]; seqs[n] = []
        elif n: seqs[n].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


by = defaultdict(list)
for h in csv.DictReader(open(hits_file), delimiter="\t"):
    by[h["sample"]].append(h)
for s, hs in by.items():
    d = out / s; d.mkdir(exist_ok=True)
    fa = d / f"{s}.windows.fa"
    if not any((d / "isescan").rglob("*.tsv")):
        seqs = read_fa(gdir / f"{s}.fa"); seen = set()
        with open(fa, "w") as f:
            for h in hs:
                st, en = int(h["region_start"]), int(h["region_end"])
                lo = max(0, st - 1 - flank); hi = min(len(seqs[h["contig"]]), en + flank)
                key = (h["contig"], lo)
                if key in seen: continue
                seen.add(key); f.write(f">{h['contig']}@{lo}\n{seqs[h['contig']][lo:hi]}\n")
        subprocess.run(["isescan.py", "--seqfile", str(fa), "--output", str(d / "isescan"), "--nthread", "1"],
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=d)
        if not any((d / "isescan").rglob("*.tsv")): print(f"[WARN] {s}: no ISEScan output", file=sys.stderr)
    print(s, file=sys.stderr)
