#!/usr/bin/env python3
"""Genomic-context analysis of blaCTX-M (default blaCTX-M-15) in E. coli assemblies.

For each gene hit: chromosome vs plasmid, relative chromosomal position, IS elements in the flanks
(ISEcp1 / IS26 credited as mediators), and other AMR genes within the window (co-integration).
Writes <outdir>/context_hits.tsv and <outdir>/summary_by_st.tsv. Standard library only.
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
                   r["Strand"], sym, r.get("Type", ""))


def read_is_hits(path, min_pident, min_len):
    hits = defaultdict(list)
    with open(path) as f:
        for line in f:
            q, s, pid, ln, qs, qe, *_ = line.rstrip("\n").split("\t")
            if float(pid) >= min_pident and int(ln) >= min_len:
                hits[q].append((min(int(qs), int(qe)), max(int(qs), int(qe)), s))
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
        isr = read_is_hits(s["is_hits"], an["is_min_pident"], an["is_min_cov_bp"])
        for contig, st, en, strand, sym, _ in amr:
            if not rx.search(sym or ""):
                continue
            lo, hi = max(1, st - win), en + win
            near_is = sorted(
                (family(n, med), "upstream" if e < st else "downstream" if b > en else "overlap",
                 max(st - e, b - en, 0))
                for b, e, n in isr.get(contig, []) if e >= lo and b <= hi)
            co = sorted({g for c, s2, e2, _, g, _ in amr
                         if c == contig and g != sym and e2 >= lo and s2 <= hi})
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
            })

    cols = list(rows[0]) if rows else []
    with open(out / "context_hits.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, cols, delimiter="\t"); w.writeheader(); w.writerows(rows)

    # Per-ST summary: gene hits, % chromosomal, mediator counts, % chromosomal hits with co-integrated AMR genes
    by = defaultdict(list)
    for r in rows: by[r["st"]].append(r)
    with open(out / "summary_by_st.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["st", "n_hits", "n_chromosomal", "pct_chromosomal", "mediators_in_chromosomal",
                    "pct_chromosomal_co_integrated"])
        for st, rs in sorted(by.items()):
            chrom = [r for r in rs if r["location"] == "chromosome"]
            m = Counter(r["mediator"] for r in chrom)
            w.writerow([st, len(rs), len(chrom), f"{100 * len(chrom) / len(rs):.1f}",
                        ";".join(f"{k}={v}" for k, v in m.most_common()) or "none",
                        f"{100 * sum(r['co_integrated'] for r in chrom) / len(chrom):.1f}" if chrom else "NA"])
    print(f"{len(rows)} hits in {len(samples)} samples -> {out}", file=sys.stderr)


if __name__ == "__main__":
    main()
