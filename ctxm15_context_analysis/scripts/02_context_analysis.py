#!/usr/bin/env python3
"""Genomic-context analysis of blaCTX-M (default blaCTX-M-15) in E. coli assemblies.

For each gene hit: chromosome vs plasmid, relative chromosomal position, IS/MGEs (MOB-recon mge.report.txt) in the flanks
(ISEcp1 / IS26 credited as mediators), and other AMR genes within the window (co-integration).
ARG barcode = sorted set of ARGs in the window; chromosomal barcodes are compared with plasmid barcodes
(identical / contained in a larger plasmid barcode / chromosome-specific) to tell isolated gene insertion
from stabilisation of a broader multidrug-resistance island.
Writes <outdir>/context_hits.tsv, barcodes.tsv and summary_by_st.tsv. Standard library only.
"""
import argparse, csv, re, sys, tomllib
from collections import Counter, defaultdict
from pathlib import Path


def contig_lengths(fasta):
    lens, name = {}, None
    for line in open(fasta):
        if line.startswith(">"):
            name = line[1:].split()[0]
            lens[name] = 0
        elif name:
            lens[name] += len(line.strip())
    return lens


def read_mobrecon(d):
    """contig_id -> (molecule_type, primary_cluster_id) from contig_report.txt."""
    out = {}
    with open(Path(d) / "contig_report.txt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r["contig_id"].split()[0]
            out[cid] = (r["molecule_type"].lower(), r.get("primary_cluster_id", "-"))
    return out


def read_amrfinder(path):
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            sym = r.get("Element symbol") or r.get("Gene symbol")
            yield (r["Contig id"].split()[0], int(r["Start"]), int(r["Stop"]),
                   r["Strand"], sym, r.get("Element type") or r.get("Type", ""))


def read_mge(path):
    """MOB-recon mge.report.txt (MOB-suite 3.x): contig -> [(start, end, mge_type)] on contig coordinates."""
    hits = defaultdict(list)
    if not Path(path).exists():
        return hits
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b, e = int(r["contig_start"]), int(r["contig_end"])
            hits[r["contig_id"].split()[0]].append((min(b, e), max(b, e), r["mge_type"]))
    return hits


def family(name, mediators):
    for m in mediators:
        if name.upper().startswith(m.upper()):
            return m
    return re.split(r"[_|]", name)[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config.toml")
    a = ap.parse_args()
    cfg = tomllib.loads(Path(a.config).read_text())
    an, inp = cfg["analysis"], cfg["input"]
    rx, win, med = re.compile(an["gene_regex"]), an["window_bp"], an["mediators"]
    out = Path(inp["outdir"]); out.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(inp["samples"]) as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    for s in samples:
        lens = contig_lengths(s["fasta"])
        mob = read_mobrecon(s["mobrecon_dir"])
        amr = list(read_amrfinder(s["amrfinder"]))
        isr = read_mge(Path(s["mobrecon_dir"]) / "mge.report.txt")
        for contig, st, en, strand, sym, _ in amr:
            if not rx.search(sym or ""):
                continue
            lo, hi = max(1, st - win), en + win
            near_is = sorted(
                (family(n, med), "upstream" if e < st else "downstream" if b > en else "overlap",
                 max(st - e, b - en, 0))
                for b, e, n in isr.get(contig, []) if e >= lo and b <= hi)
            arg_hits = [(s2, e2, g) for c, s2, e2, _, g, t in amr
                        if c == contig and g != sym and e2 >= lo and s2 <= hi and t in ("AMR", "")]
            co = sorted({g for _, _, g in arg_hits})
            # inserted region: target merged with neighbouring ARGs / IS hits separated by <= region_gap_bp
            feats = sorted([(st, en)] + [(a, b) for a, b, _ in arg_hits] +
                           [(b, e) for b, e, _ in isr.get(contig, []) if e >= lo and b <= hi])
            rs, re_ = st, en
            for _ in feats:
                for b, e in feats:
                    if b <= re_ + an["region_gap_bp"] and e >= rs - an["region_gap_bp"]:
                        rs, re_ = min(rs, b), max(re_, e)
            mols = mob.get(contig, ("unknown", "-"))
            clen = lens.get(contig, 0)
            nearest = {m: min((d for f_, _, d in near_is if f_ == m), default=None) for m in med}
            mediator = next((m for m in med if nearest[m] is not None), "none")
            rows.append({
                "sample": s["sample"], "st": s.get("st", "NA"), "gene": sym, "contig": contig,
                "contig_len": clen, "start": st, "end": en, "strand": strand,
                "location": mols[0], "plasmid_cluster": mols[1],
                "rel_position": f"{(st + en) / 2 / clen:.4f}" if clen else "NA",
                "mediator": mediator,
                "is_in_window": ";".join(f"{f_}:{side}:{d}" for f_, side, d in near_is) or "none",
                "co_resistance_genes": ";".join(co) or "none",
                "co_integrated": int(bool(co)),
                "region_start": rs, "region_end": re_,
                "barcode": ",".join(sorted({sym, *co})),
                "barcode_class": "plasmid" if mols[0] == "plasmid" else "",
            })

    # Barcode classes: chromosomal barcode vs plasmid barcodes across all genomes
    pl = {frozenset(r["barcode"].split(",")) for r in rows if r["location"] == "plasmid"}
    for r in rows:
        if r["location"] == "plasmid":
            continue
        b = frozenset(r["barcode"].split(","))
        r["barcode_class"] = ("identical_to_plasmid" if b in pl else
                              "contained_in_larger_plasmid" if any(b < p for p in pl) else
                              "chromosome_specific") if r["location"] == "chromosome" else r["location"]

    cols = list(rows[0]) if rows else []
    with open(out / "context_hits.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

    # Barcodes by location, with recurrence across genomes, variants and STs
    bc = defaultdict(list)
    for r in rows: bc[(r["location"], r["barcode"])].append(r)
    with open(out / "barcodes.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["location", "barcode", "n_args", "n_genomes", "variants", "sts", "barcode_class"])
        for (loc, b), rs in sorted(bc.items(), key=lambda kv: -len({r["sample"] for r in kv[1]})):
            w.writerow([loc, b, len(b.split(",")), len({r["sample"] for r in rs}),
                        ";".join(sorted({r["gene"] for r in rs})), ";".join(sorted({r["st"] for r in rs})),
                        rs[0]["barcode_class"]])

    # Per-ST summary, gene level (each copy counts) and genome level (>=1 chromosomal copy): gene hits, % chromosomal, mediator counts, % chromosomal hits with co-integrated AMR genes
    by = defaultdict(list)
    for r in rows: by[r["st"]].append(r)
    with open(out / "summary_by_st.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["st", "n_hits", "n_chromosomal", "pct_chromosomal", "mediators_in_chromosomal",
                    "pct_chromosomal_co_integrated", "n_genomes", "n_genomes_chromosomal"])
        for st, rs in sorted(by.items()):
            chrom = [r for r in rs if r["location"] == "chromosome"]
            m = Counter(r["mediator"] for r in chrom)
            w.writerow([st, len(rs), len(chrom), f"{100 * len(chrom) / len(rs):.1f}",
                        ";".join(f"{k}={v}" for k, v in m.most_common()) or "none",
                        f"{100 * sum(r['co_integrated'] for r in chrom) / len(chrom):.1f}" if chrom else "NA",
                        len({r["sample"] for r in rs}), len({r["sample"] for r in chrom})])
    print(f"{len(rows)} hits in {len(samples)} samples -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
