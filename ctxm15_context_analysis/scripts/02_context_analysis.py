#!/usr/bin/env python3
"""Genomic-context analysis of blaCTX-M (default blaCTX-M-15) in E. coli assemblies.

Tools (config [analysis]): arg_tool = "resfinder" (ResFinder 4 tab output, default) | "amrfinder" (AMRFinderPlus);
location_tool = "platon" (Platon + size rule: contigs > plasmid_size_rule are chromosomal, default) | "mobrecon" (MOB-recon);
is_tool = "isfinder" (ABRicate vs ISfinder, default) | "mobrecon" (MOB-recon mge.report; no IS26 entry, fewer ISEcp1 hits).
If MOB-recon output is also given it supplies plasmid clusters and IS/MGEs, and any chromosome/plasmid disagreement is flagged.
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
    """contig_id -> (molecule_type, primary_cluster_id, circular) from contig_report.txt."""
    out = {}
    with open(Path(d) / "contig_report.txt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r["contig_id"].split()[0]
            out[cid] = (r["molecule_type"].lower(), r.get("primary_cluster_id", "-"),
                        r.get("circularity_status", "").lower().startswith("circular"))
    return out


def read_amrfinder(path):
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("Element subtype") == "POINT":      # --organism point mutations are not ARG context
                continue
            sym = r.get("Element symbol") or r.get("Gene symbol")
            yield (r["Contig id"].split()[0], int(r["Start"]), int(r["Stop"]),
                   r["Strand"], sym, r.get("Element type") or r.get("Type", ""))


TIES = {}      # (contig, start, end) -> alleles that tied for that locus (ResFinder lists every equally good allele)


def read_resfinder(path, min_id, min_cov):
    """ResFinder 4 ResFinder_results_tab.txt -> same tuples as read_amrfinder (strand unknown: 'NA').
    ResFinder lists every allele that ties at a locus; overlapping hits (> 50% of the shorter) are collapsed to the best
    identity x coverage hit, first listed on ties, and the tied alleles are kept in TIES."""
    hits = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            idt, cov = float(r["Identity"]), float(r["Coverage"])
            if idt < min_id or cov < min_cov:
                continue
            a, b = map(int, r["Position in contig"].split(".."))
            hits.append((r["Contig"].split()[0], min(a, b), max(a, b), r["Resistance gene"], idt * cov))
    kept = []
    for h in sorted(hits, key=lambda h: -h[4]):                 # stable: file order breaks ties
        for k in kept:
            ov = min(h[2], k[2]) - max(h[1], k[1])
            if k[0] == h[0] and ov > 0.5 * min(h[2] - h[1], k[2] - k[1]):
                if h[4] == k[4] and h[3] != k[3]: TIES.setdefault((k[0], k[1], k[2]), {k[3]}).add(h[3])
                break
        else:
            kept.append(h)
    for c, st, en, g, _ in sorted(kept, key=lambda h: (h[0], h[1])):
        yield (c, st, en, "NA", g, "AMR")


def read_platon(d):
    """Platon output dir -> (plasmid contig ids, circular plasmid contig ids), or None if not run.
    Platon writes no table when no contig passes its filters (e.g. one complete chromosome): that means no plasmid contigs."""
    d = Path(d)
    tsv = next(iter(sorted(d.glob("*.tsv"))), None)
    if tsv is None:
        log = next(iter(sorted(d.glob("*.log"))), None)
        return (set(), set()) if log and "no potential plasmid contigs" in log.read_text() else None
    ids, circ = set(), set()
    with open(tsv) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            cid = r["ID"].split()[0]
            ids.add(cid)
            if r.get("Circular", "").lower() == "yes": circ.add(cid)
    return ids, circ


def gbff_circular(path):
    """Contig ids (LOCUS name and VERSION) of circular records in a GBFF."""
    out, cur = set(), None
    if not path or not Path(path).exists():
        return out
    for line in open(path):
        if line.startswith("LOCUS"):
            cur = line.split()[1]
            if "circular" in line: out.add(cur)
        elif line.startswith("VERSION") and cur in out:
            out.add(line.split()[1])
        elif line.startswith("FEATURES"):
            cur = None
    return out


def read_mge(path):
    """MOB-recon mge.report.txt (MOB-suite 3.x): contig -> [(start, end, mge_type)] on contig coordinates."""
    hits = defaultdict(list)
    if not Path(path).exists():
        return hits
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b, e = int(r["contig_start"]), int(r["contig_end"])
            name = r["mge_type"]
            if r.get("mge_subtype") == "IS6" and not name.startswith("IS26"):
                name = "IS26"            # MOB-suite's MGE set has no IS26 entry; its IS6-family copies (IS6100, ...) count as IS26-family
            hits[r["contig_id"].split()[0]].append((min(b, e), max(b, e), name))
    return hits


ISFINDER_ALIAS = {"ISEc9": "ISEcp1"}   # ISfinder lists ISEc9 as a 100%-identical copy of ISEcp1; ABRicate prints one name per locus


def read_abricate_is(path):
    """ABRicate output against an ISfinder database -> contig -> [(start, end, name)] (ISEc9 reported as ISEcp1)."""
    hits = defaultdict(list)
    if not path or not Path(path).exists():
        return hits
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            n = r["GENE"].split("_")[0]
            b, e = int(r["START"]), int(r["END"])
            hits[r["SEQUENCE"].split()[0]].append((min(b, e), max(b, e), ISFINDER_ALIAS.get(n, n)))
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

    arg_tool, loc_tool = an.get("arg_tool", "amrfinder"), an.get("location_tool", "mobrecon")
    size_rule = an.get("plasmid_size_rule", 500_000)
    is_tool = an.get("is_tool", "mobrecon")
    rows = []
    with open(inp["samples"]) as f:
        samples = list(csv.DictReader(f, delimiter="\t"))
    for s in samples:
        lens = contig_lengths(s["fasta"])
        mob = read_mobrecon(s["mobrecon_dir"]) if s.get("mobrecon_dir") else {}
        if arg_tool == "resfinder":
            amr = list(read_resfinder(s["resfinder"], an.get("resfinder_min_identity", 80.0), an.get("resfinder_min_coverage", 80.0)))
        else:
            amr = list(read_amrfinder(s["amrfinder"]))
        if is_tool == "isfinder":
            isr = read_abricate_is(s.get("isfinder"))
        else:
            isr = read_mge(Path(s["mobrecon_dir"]) / "mge.report.txt") if s.get("mobrecon_dir") else defaultdict(list)
        plat = read_platon(s["platon"]) if loc_tool == "platon" else None
        if loc_tool == "platon" and plat is None:
            print(f"[WARN] {s['sample']}: no Platon output, locations unknown", file=sys.stderr)
        circular = gbff_circular(s.get("gbff")) | (plat[1] if plat else set()) | {c for c, v in mob.items() if v[2]}
        for contig, st, en, strand, sym, _ in amr:
            if not rx.search(sym or ""):
                continue
            lo, hi = max(1, st - win), en + win
            near_is = sorted(
                (family(n, med), "upstream" if e < st else "downstream" if b > en else "overlap",
                 max(st - e, b - en, 0))
                for b, e, n in isr.get(contig, []) if e >= lo and b <= hi)
            # inserted region: the target plus everything within window_bp, extended by ARGs / IS chained on with gaps <= region_gap_bp (capped at region_max_bp each side).
            # The ARG barcode is the set of ARGs inside this region, so islands longer than the flank window are not cut off.
            span = an.get("region_max_bp", 50_000)
            feats = [(s2, e2, g, True) for c, s2, e2, _, g, t in amr
                     if c == contig and (g != sym or (s2, e2) != (st, en)) and t in ("AMR", "")]
            feats += [(b, e, n, False) for b, e, n in isr.get(contig, [])]
            feats = [f for f in feats if f[1] >= st - span and f[0] <= en + span]
            rs, re_ = st, en
            for b, e, _, _ in feats:                       # everything inside the flank window belongs to the region...
                if e >= lo and b <= hi: rs, re_ = min(rs, b), max(re_, e)
            changed = True
            while changed:                                 # ...and so does anything chained on beyond it
                changed = False
                for b, e, _, _ in feats:
                    if b <= re_ + an["region_gap_bp"] and e >= rs - an["region_gap_bp"] and (b < rs or e > re_):
                        rs, re_ = min(rs, b), max(re_, e); changed = True
            arg_hits = [(b, e, g) for b, e, g, is_arg in feats if is_arg and e >= rs and b <= re_]
            co = sorted({g for _, _, g in arg_hits if g != sym})
            clen = lens.get(contig, 0)
            mobm = mob.get(contig)
            if loc_tool == "platon":
                if plat is None: loc = "unknown"
                elif clen > size_rule: loc = "chromosome"            # paper rule: large contigs are chromosomal
                else: loc = "plasmid" if contig in plat[0] else "chromosome"
                mols = (loc, mobm[1] if mobm else "-", contig in circular)
            else:
                mols = mobm or ("unknown", "-", False)
                mols = (mols[0], mols[1], contig in circular)
            disagree = ("mobrecon:" + mobm[0]) if (loc_tool == "platon" and mobm and mobm[0] != mols[0]) else "-"
            # unreliable-context flags: only linear contigs can be cut short by their ends
            q = []
            if not mols[2]:
                if clen < an["min_contig_bp"]: q.append("short_contig")
                if lo <= 1 or hi >= clen: q.append("truncated_flank")
            nearest = {m: min((d for f_, _, d in near_is if f_ == m), default=None) for m in med}
            mediator = next((m for m in med if nearest[m] is not None), "none")
            rows.append({
                "sample": s["sample"], "st": s.get("st", "NA"), "gene": sym, "contig": contig,
                "contig_len": clen, "start": st, "end": en, "strand": strand,
                "allele_ties": ";".join(sorted(TIES.get((contig, st, en), []))) or "-",
                "location": mols[0], "plasmid_cluster": mols[1], "location_disagreement": disagree,
                "rel_position": f"{(st + en) / 2 / clen:.4f}" if clen else "NA",
                "mediator": mediator,
                "is_in_window": ";".join(f"{f_}:{side}:{d}" for f_, side, d in near_is) or "none",
                "co_resistance_genes": ";".join(co) or "none",
                "co_integrated": int(bool(co)),
                "region_start": rs, "region_end": re_,
                "context_quality": ",".join(q) or "ok",
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
