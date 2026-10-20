#!/usr/bin/env python3
"""Collect per-genome typing (mlst, ClermonTyping, ECTyper, FimTyper) and assign ST131 subclades.

Reads samples.tsv columns mlst, clermont, ectyper, fimtyper (any may be blank). Writes results/typing.tsv.
Parsing is header-tolerant because output layouts differ between tool versions; check typing.tsv.
ST131 subclade (heuristic from fimH / O-type / blaCTX-M-15; replace with the paper's reference 47 definitions
once available): fimH41 -> A, fimH22 -> B, fimH30 -> C2 if the genome carries blaCTX-M-15 else C1, any other fimH -> Other.
Needs results/context_hits.tsv for the blaCTX-M-15 check (run 02 with gene_regex "^blaCTX-M" for all variants).
"""
import csv, re, sys, tomllib
from pathlib import Path


def rows(path):
    return [] if not path or not Path(path).exists() else list(csv.reader(open(path), delimiter="\t"))


def get_mlst(path):                       # file, scheme, ST, alleles...
    r = rows(path)
    return f"ST{r[0][2]}" if r and len(r[0]) > 2 and r[0][2] not in ("-", "") else "ND"


def get_clermont(path):                   # ClermonTyping *_phylogroups.txt: no header; col 5 = phylogroup (col 6 = Mash group)
    r = [x for x in rows(path) if x]
    if not r: return "ND"
    hdr = next((k for k, h in enumerate(r[0]) if "phylogroup" in h.lower()), None)
    if hdr is not None and len(r) > 1: return r[-1][hdr] or "ND"
    return (r[-1][4] if len(r[-1]) > 4 else r[-1][-1]) or "ND"


def get_ectyper(path):                    # columns O-type, H-type
    r = rows(path)
    if len(r) < 2: return "ND", "ND"
    h = {k: i for i, k in enumerate(r[0])}
    return (r[1][h["O-type"]] if "O-type" in h else "ND", r[1][h["H-type"]] if "H-type" in h else "ND")


def get_fimh(path):                       # FimTyper results_tab.txt: "FimH type: FimH30" (also accepts a Fimtype table)
    if not path or not Path(path).exists(): return "ND"
    m = re.search(r"fimh\s*(\d+)", Path(path).read_text(), re.I)
    return f"fimH{m.group(1)}" if m else "ND"


def subclade(st, fimh, ctxm15):
    if st != "ST131": return "NA"
    if fimh == "ND": return "ND"
    return {"fimH41": "A", "fimH22": "B"}.get(fimh, ("C2" if ctxm15 else "C1") if fimh == "fimH30" else "Other")


def main():
    cfg = tomllib.loads(Path(sys.argv[1] if len(sys.argv) > 1 else "config.toml").read_text())
    out = Path(cfg["input"]["outdir"]); out.mkdir(parents=True, exist_ok=True)
    c15 = set()
    hp = out / "context_hits.tsv"
    if hp.exists():
        c15 = {r["sample"] for r in csv.DictReader(open(hp), delimiter="\t") if r["gene"] == "blaCTX-M-15"}
    with open(out / "typing.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["sample", "st", "phylogroup", "o_type", "h_type", "fimh", "carries_ctxm15", "st131_subclade"])
        for s in csv.DictReader(open(cfg["input"]["samples"]), delimiter="\t"):
            st, fh = get_mlst(s.get("mlst")), get_fimh(s.get("fimtyper"))
            o, h = get_ectyper(s.get("ectyper"))
            w.writerow([s["sample"], st, get_clermont(s.get("clermont")), o, h, fh,
                        int(s["sample"] in c15), subclade(st, fh, s["sample"] in c15)])


if __name__ == "__main__":
    main()
