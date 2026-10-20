#!/usr/bin/env python3
"""Extract blaCTX-M windows from the GBFF files and align them per variant with Clinker.

For each variant, windows = inserted region (context_hits.tsv region_start/end) +- flank_bp, cut from the
sample GBFF, reverse-complemented so blaCTX-M is on the forward strand, written to
<outdir>/clinker/<variant>/<sample>.gbk. Hits are ordered by barcode recurrence (most frequent structures
first) and capped at max_per_variant. With --run, `clinker <variant>/*.gbk -p <variant>.html -o <variant>.csv`
is executed for each variant. Needs Biopython (and clinker for --run).
"""
import argparse, csv, re, subprocess, sys, tomllib
from collections import Counter, defaultdict
from pathlib import Path

from Bio import SeqIO


def gene_strand(rec, start, end):
    """Strand of the CDS that overlaps the gene most (ResFinder hits carry no strand); '+' if none.
    Ties go to the shorter CDS: some GBFFs carry a CDS spanning a whole contig that would otherwise win."""
    best, strand = (0, 0), "+"
    for f in rec.features:
        if f.type == "CDS":
            ov = min(int(f.location.end), end) - max(int(f.location.start) + 1, start) + 1
            key = (ov, -(int(f.location.end) - int(f.location.start)))
            if ov > 0 and key > best: best, strand = key, "-" if f.location.strand == -1 else "+"
    return strand


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml")
    ap.add_argument("--run", action="store_true", help="run clinker on each variant directory")
    a = ap.parse_args()
    cfg = tomllib.loads(Path(a.config).read_text())
    cc, inp = cfg["clinker"], cfg["input"]
    out = Path(inp["outdir"]); root = out / "clinker"

    gbff = {r["sample"]: r.get("gbff") for r in csv.DictReader(open(inp["samples"]), delimiter="\t")}
    hits = [h for h in csv.DictReader(open(out / "context_hits.tsv"), delimiter="\t")
            if cc["location"] in ("all", h["location"]) and gbff.get(h["sample"])]
    freq = Counter((h["gene"], h["barcode"]) for h in hits)
    by = defaultdict(list)
    for h in hits: by[h["gene"]].append(h)

    recs = {}
    for variant, hs in by.items():
        d = root / variant; d.mkdir(parents=True, exist_ok=True)
        hs.sort(key=lambda h: (-freq[(h["gene"], h["barcode"])], h["sample"]))
        seen, rows = set(), []
        for h in hs:
            if len(rows) >= cc["max_per_variant"]:
                break
            if h["sample"] in seen:                      # one window per genome so Clinker rows are genomes
                continue
            seen.add(h["sample"])
            if h["sample"] not in recs:
                recs[h["sample"]] = {r.id: r for r in SeqIO.parse(gbff[h["sample"]], "genbank")}
            rec = recs[h["sample"]].get(h["contig"])
            if rec is None:
                print(f"[WARN] {h['sample']}: contig {h['contig']} not in GBFF", file=sys.stderr); continue
            s = max(0, int(h["region_start"]) - 1 - cc["flank_bp"])
            e = min(len(rec), int(h["region_end"]) + cc["flank_bp"])
            sub = rec[s:e]
            strand = h["strand"] if h["strand"] in ("+", "-") else gene_strand(rec, int(h["start"]), int(h["end"]))
            if strand == "-":
                sub = sub.reverse_complement(id=True, name=True, description=True, annotations=True)
            name = re.sub(r"\W", "_", h["sample"])
            sub.id = sub.name = name[:16]
            sub.description = f"{h['sample']} {h['contig']}:{s + 1}-{e} {h['location']} barcode={h['barcode']}"
            sub.annotations["molecule_type"] = "DNA"
            SeqIO.write(sub, d / f"{name}.gbk", "genbank")
            rows.append([h["sample"], h["contig"], s + 1, e, strand, h["mediator"], h["barcode"]])
        with open(d / "manifest.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["sample", "contig", "window_start", "window_end", "strand", "mediator", "barcode"]); w.writerows(rows)
        print(f"{variant}: {len(rows)} windows -> {d}", file=sys.stderr)
        if a.run and len(rows) >= 2:
            subprocess.run(["clinker", *map(str, sorted(d.glob("*.gbk"))), "-p", str(root / f"{variant}.html"),
                            "-o", str(root / f"{variant}.csv"), "-f"], check=True)


if __name__ == "__main__":
    main()
