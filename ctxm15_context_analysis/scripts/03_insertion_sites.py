#!/usr/bin/env python3
"""Map chromosomal blaCTX-M insertion sites onto the E. coli K-12 0-100 min scale.

Input: results/context_hits.tsv (chromosomal rows) + per-sample GBFF (column `gbff` of samples.tsv; LOCUS or VERSION ids must match contig ids).
Flanking CDS on each side of the inserted region are the anchors. Mobile-element genes (transposase, integrase, recombinase,
IS) are skipped, circular chromosomes wrap across the origin, and each anchor is mapped to K-12 by
  mode = "homology" (default): BLASTp of the anchor protein against the K-12 proteome (best hit >= min_identity and
                   >= min_coverage, unique: best bitscore >= ambiguity_ratio x second best), falling back to the gene symbol / locus tag;
  mode = "name":   gene symbol / locus tag / product only (the first version; fails for genes with no symbol).
Scanning goes outward (up to max_scan CDS per flank) until n_anchors anchors are mapped; sites whose nearest anchor lies more than island_skip CDS away
are flagged site_quality = island_flank (genomic islands absent from K-12 push the anchors outwards). The longest contiguous anchor chain per
flank (consecutive gap <= max_chain_gap_min, >= min_chain_genes genes) fixes the position:
  both flanks valid and junction gap <= max_junction_gap_min -> circular midpoint of the nearest anchors
  junction gap larger -> flank with the longer chain (equal chains -> ND)
  one flank valid -> nearest anchor of that flank;  none -> ND
Writes insertion_sites.tsv and site_histogram.tsv (genomes per 0.1 min). Needs blastp/makeblastdb on PATH for the homology mode.
Usage: 03_insertion_sites.py [config.toml] [--mode homology|name] [--out-suffix _name]
"""
import argparse, csv, re, subprocess, sys, tempfile, tomllib
from collections import Counter
from pathlib import Path

MOBILE = re.compile(r"transposase|integrase|recombinase|insertion (element|sequence)|\bIS\d|\bIS[A-Z][a-z]|mobile element|resolvase", re.I)


def parse_gbff(gbff):
    """Return ({record: [CDS dict]}, {circular record ids}); CDS dict = start, end, strand, keys, product, prot."""
    recs, circular, rec = {}, set(), None
    feat, in_feat, q = None, False, None

    def close():
        nonlocal feat, q
        if feat and feat["key"] == "CDS":
            nums = list(map(int, re.findall(r"\d+", feat["loc"])))
            if nums:
                quals = {k: "".join(v).strip('"') for k, v in feat["q"].items()}
                keys = {quals[k].strip().lower() for k in ("gene", "locus_tag", "product") if quals.get(k)}
                recs[rec].append({"start": min(nums), "end": max(nums), "strand": -1 if feat["loc"].startswith("complement") else 1,
                                  "keys": keys, "product": quals.get("product", ""), "prot": quals.get("translation", "").replace(" ", ""),
                                  "id": quals.get("locus_tag") or quals.get("protein_id") or f"cds{len(recs[rec])}"})
        feat, q = None, None

    for line in open(gbff):
        if line.startswith("LOCUS"):
            close(); in_feat = False
            rec = line.split()[1]; recs[rec] = []
            if "circular" in line: circular.add(rec)
        elif line.startswith("VERSION") and rec:
            recs[line.split()[1]] = recs[rec]
            if rec in circular: circular.add(line.split()[1])
        elif line.startswith("FEATURES"):
            in_feat = True
        elif line.startswith("ORIGIN") or line.startswith("CONTIG"):
            close(); in_feat = False
        elif in_feat:
            if line[5:6].strip():                                  # new feature: key at col 6, location at col 22
                close(); feat = {"key": line[5:21].strip(), "loc": line[21:].strip(), "q": {}}
            elif feat:
                t = line[21:].rstrip("\n")
                if t.lstrip().startswith("/"):
                    m = re.match(r"\s*/([^=]+)(?:=(.*))?", t)
                    q = m.group(1); feat["q"][q] = [(m.group(2) or "").strip()]
                elif q is None:
                    feat["loc"] += t.strip()                       # wrapped location (join/complement)
                else:
                    feat["q"][q].append(t.strip())
    close()
    return recs, circular


def circ(a, b):                                                    # signed circular difference b-a on a 100 scale
    return (b - a + 50) % 100 - 50


def chain(minutes, max_gap, min_genes):
    """Longest run of consecutive anchors (nearest first) with gaps <= max_gap; ties -> nearest."""
    best, run = [], []
    for m in minutes:
        run = run + [m] if run and abs(circ(run[-1], m)) <= max_gap else [m]
        if len(run) > len(best): best = run
    return best if len(best) >= min_genes else []


def place(mu, md, sc):
    """(position in min or None, method, upstream chain, downstream chain) from the mapped anchor minutes of each flank."""
    cu = chain(mu, sc["max_chain_gap_min"], sc["min_chain_genes"])
    cd = chain(md, sc["max_chain_gap_min"], sc["min_chain_genes"])
    pos, how = None, "ND"
    if cu and cd and abs(circ(cu[0], cd[0])) <= sc["max_junction_gap_min"]:
        pos, how = (cu[0] + circ(cu[0], cd[0]) / 2) % 100, "both_flanks"
    elif cu and cd:
        if len(cu) != len(cd):
            pos, how = (cu if len(cu) > len(cd) else cd)[0], "longer_flank"
    elif cu or cd:
        pos, how = (cu or cd)[0], "single_flank"
    return pos, how, cu, cd


def flank_candidates(cds, rs, re_, circular, max_scan):
    """(upstream, downstream) CDS lists ordered nearest first, outside [rs, re_]; wrap across the origin if circular."""
    cl = sorted(cds, key=lambda c: c["start"])
    n = len(cl)
    up_i = [i for i, c in enumerate(cl) if c["end"] < rs]
    dn_i = [i for i, c in enumerate(cl) if c["start"] > re_]
    up = [cl[i] for i in reversed(up_i)]
    dn = [cl[i] for i in dn_i]
    if circular and n:
        up += [cl[i] for i in reversed(dn_i)]                      # continue past the start of the record: the end comes next
        dn += [cl[i] for i in up_i]
    return up[:max_scan], dn[:max_scan]


def build_k12(k12_gbff):
    recs, _ = parse_gbff(k12_gbff)
    cds = list({id(c): c for cs in recs.values() for c in cs}.values())     # records are indexed under LOCUS and VERSION: count each CDS once
    glen = max(c["end"] for c in cds)
    byname, byid = {}, {}
    for c in cds:
        mid = (c["start"] + c["end"]) / 2 / glen * 100
        byid[c["id"]] = mid
        for k in c["keys"]: byname.setdefault(k, mid)
    return cds, glen, byname, byid


def blast_map(queries, k12_cds, byid, sc, tmp):
    """{query id: K-12 minute or None} by unique best BLASTp hit."""
    db = Path(tmp) / "k12.faa"
    with open(db, "w") as f:
        for c in k12_cds:
            if c["prot"]: f.write(f">{c['id']}\n{c['prot']}\n")
    subprocess.run(["makeblastdb", "-in", str(db), "-dbtype", "prot", "-out", str(db)], check=True, capture_output=True)
    qf = Path(tmp) / "q.faa"
    with open(qf, "w") as f:
        for qid, prot in queries.items(): f.write(f">{qid}\n{prot}\n")
    out = subprocess.run(["blastp", "-query", str(qf), "-db", str(db), "-evalue", "1e-10", "-max_target_seqs", "5",
                          "-outfmt", "6 qseqid sseqid pident qcovs bitscore"], check=True, capture_output=True, text=True).stdout
    hits = {}
    for line in out.splitlines():
        q, s, pid, qcov, bits = line.split("\t")
        hits.setdefault(q, []).append((float(bits), s, float(pid), float(qcov)))
    res = {}
    for q in queries:
        hs = sorted(hits.get(q, []), reverse=True)
        if hs and hs[0][2] >= sc["min_identity"] and hs[0][3] >= sc["min_coverage"] and \
                (len(hs) == 1 or hs[0][0] >= sc["ambiguity_ratio"] * hs[1][0]):
            res[q] = byid[hs[0][1]]
        else:
            res[q] = None
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config", nargs="?", default="config.toml")
    ap.add_argument("--mode", choices=["homology", "name"], default=None)
    ap.add_argument("--out-suffix", default="")
    a = ap.parse_args()
    cfg = tomllib.loads(Path(a.config).read_text())
    sc, inp = cfg["sites"], cfg["input"]
    mode = a.mode or sc.get("mode", "homology")
    sc.setdefault("max_scan", 80); sc.setdefault("island_skip", 10); sc.setdefault("min_identity", 80.0); sc.setdefault("min_coverage", 80.0); sc.setdefault("ambiguity_ratio", 1.1)
    out = Path(inp["outdir"])

    k12_cds, glen, kname, kid = build_k12(sc["k12_gbff"])
    gbff = {r["sample"]: r.get("gbff") for r in csv.DictReader(open(inp["samples"]), delimiter="\t")}
    hits = [h for h in csv.DictReader(open(out / "context_hits.tsv"), delimiter="\t")
            if h["location"] == "chromosome" and gbff.get(h["sample"])]

    rows, parsed, mapped = [], {}, {}
    for h in hits:
        s = h["sample"]
        if s not in parsed:
            parsed[s] = parse_gbff(gbff[s])
        recs, circ_ids = parsed[s]
        cds = recs.get(h["contig"], [])
        rs, re_ = int(h["region_start"]), int(h["region_end"])
        up, dn = flank_candidates(cds, rs, re_, h["contig"] in circ_ids, sc["max_scan"])
        h["_cands"] = (up, dn)

    if mode == "homology":
        with tempfile.TemporaryDirectory() as tmp:
            queries = {}
            for h in hits:
                for side in h["_cands"]:
                    for c in side:
                        if c["prot"] and not MOBILE.search(c["product"]): queries[f"{h['sample']}|{c['id']}"] = c["prot"]
            mapped = blast_map(queries, k12_cds, kid, sc, tmp)

    def anchors(side, sample):
        """(minutes of the first n_anchors mappable non-mobile CDS, nearest first; CDS skipped before the first one)."""
        ms, skipped = [], None
        for i, c in enumerate(side):
            if MOBILE.search(c["product"]): continue
            m = mapped.get(f"{sample}|{c['id']}") if mode == "homology" else None
            if m is None:                                          # gene symbol / locus tag / product
                m = next((kname[k] for k in c["keys"] if k in kname), None)
            if m is not None:
                ms.append(m)
                if skipped is None: skipped = i
            if len(ms) >= sc["n_anchors"]: break
        return ms, skipped

    for h in hits:
        up, dn = h["_cands"]
        mu, su = anchors(up, h["sample"]); md, sd = anchors(dn, h["sample"])
        pos, how, cu, cd = place(mu, md, sc)
        rows.append({"sample": h["sample"], "st": h["st"], "gene": h["gene"],
                     "position_min": "ND" if pos is None else f"{round(pos, 1) % 100:.1f}",
                     "method": how, "n_up_chain": len(cu), "n_down_chain": len(cd),
                     "up_anchor_min": f"{cu[0]:.2f}" if cu else "NA",
                     "down_anchor_min": f"{cd[0]:.2f}" if cd else "NA",
                     "up_skipped": "NA" if su is None else su, "down_skipped": "NA" if sd is None else sd,
                     # nearest mapped gene far from the region = the insertion lies in / next to a genomic island absent from K-12;
                     # the position is then only accurate to about the island size (~1 min per 46 kb)
                     "site_quality": "island_flank" if any(x is not None and x > sc["island_skip"] for x in (su, sd)) else "ok",
                     "mediator": h["mediator"], "anchor_mode": mode})

    sfx = a.out_suffix
    with open(out / f"insertion_sites{sfx}.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]) if rows else [], delimiter="\t"); w.writeheader(); w.writerows(rows)
    hist = Counter(r["position_min"] for r in rows if r["position_min"] != "ND")
    with open(out / f"site_histogram{sfx}.tsv", "w") as f:
        f.write("position_min\tn_genomes\n")
        for p, n in sorted(hist.items(), key=lambda kv: float(kv[0])): f.write(f"{p}\t{n}\n")
    c = Counter(r["method"] for r in rows)
    print(f"[{mode}] {len(rows)} chromosomal insertions: " + ", ".join(f"{k} {v}" for k, v in c.most_common()), file=sys.stderr)


if __name__ == "__main__":
    main()
