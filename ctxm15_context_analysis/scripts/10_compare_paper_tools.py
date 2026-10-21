#!/usr/bin/env python3
"""Compare three views of blaCTX-M location, mediator and typing on the same genomes:
  paper table   = Zenodo summary table (the authors' calls, ground truth for copy numbers / ST131 clade)
  ours          = AMRFinderPlus + MOB-recon + ECTyper/FimTyper pipeline (steps 01-05)
  paper tools   = ABRicate (ResFinder / ISfinder / FimTyper / EcOH) + Platon with the 500 kb size rule
Usage: 10_compare_paper_tools.py <config.toml> <subset.tsv> <run_paper_dir> [<fasttree_newick>]
Writes <run_paper_dir>/comparison/{ctxm_hits_paper_tools.tsv, per_genome_comparison.tsv, tree_rf.tsv, summary.txt}.
"""
import csv, re, sys, tomllib
from collections import defaultdict
from pathlib import Path

WIN = 10000                      # same flank window as step 02
# ISfinder lists ISEc9 as a 100%-identical duplicate of ISEcp1 and ABRicate reports only one name per locus (ISEc9 here)
MEDIATORS = {"ISEcp1": ("ISEcp1", "ISEc9"), "IS26": ("IS26",)}
SIZE_RULE = 500_000              # paper: contigs > 500 kb are chromosomal / reviewed manually


def contig_lengths(fa):
    L, n = {}, None
    for line in open(fa):
        if line.startswith(">"): n = line[1:].split()[0]; L[n] = 0
        elif n: L[n] += len(line.strip())
    return L


def abricate(path):
    if not Path(path).exists(): return []
    return list(csv.DictReader(open(path), delimiter="\t"))


def platon_plasmid_contigs(d, sample):
    p = Path(d) / "platon" / sample / f"{sample}.tsv"
    if p.exists(): return {r["ID"].split()[0] for r in csv.DictReader(open(p), delimiter="\t")}
    log = p.with_suffix(".log")      # Platon writes no table when no contig passes its filters (e.g. a single complete chromosome)
    if log.exists() and "no potential plasmid contigs" in log.read_text(): return set()
    return None


def collapse(hits):
    """Keep one blaCTX-M hit per locus (ResFinder lists several overlapping alleles): best identity*coverage."""
    out = []
    for h in sorted(hits, key=lambda h: -(float(h["%IDENTITY"]) * float(h["%COVERAGE"]))):
        s, e = int(h["START"]), int(h["END"])
        if not any(o["SEQUENCE"] == h["SEQUENCE"] and min(e, int(o["END"])) - max(s, int(o["START"])) > 0.5 * (e - s) for o in out):
            out.append(h)
    return out


def variant(g):
    m = re.match(r"(blaCTX-M-\d+)", g)
    return m.group(1) if m else g.split("_")[0]


def paper_counts(v):
    return {m.group(1): int(m.group(2)) for m in re.finditer(r"(blaCTX-M-[\w]+) \(n=(\d+)\)", v)}


def bipartitions(tree):
    tips = {t.name for t in tree.get_terminals()}
    out = set()
    for cl in tree.get_nonterminals():
        s = frozenset(t.name for t in cl.get_terminals())
        if 1 < len(s) < len(tips) - 1: out.add(min(s, frozenset(tips - s), key=lambda x: sorted(x)))
    return out


def main():
    cfg = tomllib.loads(Path(sys.argv[1]).read_text()); sub = list(csv.DictReader(open(sys.argv[2]), delimiter="\t"))
    rp = Path(sys.argv[3]); out = rp / "comparison"; out.mkdir(exist_ok=True)
    res = Path(cfg["input"]["outdir"])
    ours = defaultdict(list)
    for h in csv.DictReader(open(res / "context_hits.tsv"), delimiter="\t"): ours[h["sample"]].append(h)
    typ = {r["sample"]: r for r in csv.DictReader(open(res / "typing.tsv"), delimiter="\t")}

    rows, hit_rows = [], []
    for r in sub:
        s = r["Strain"]; L = contig_lengths(f"test_data/genomes/{s}.fa")
        pl = platon_plasmid_contigs(rp, s)
        ctx = collapse([h for h in abricate(rp / "abricate" / f"{s}.resfinder.tsv") if h["GENE"].startswith("blaCTX-M")])
        isf = abricate(rp / "abricate" / f"{s}.isfinder.tsv")
        pt = {"chromosome": 0, "plasmid": 0}; pv_chrom = set(); med = set(); review = 0
        for h in ctx:
            c, st, en = h["SEQUENCE"], int(h["START"]), int(h["END"])
            if pl is None: loc = "NA"
            elif L.get(c, 0) > SIZE_RULE: loc = "chromosome"; review += c in pl      # large contig called plasmid by Platon: flagged
            elif c in pl: loc = "plasmid"
            else: loc = "chromosome"
            near = [(i["GENE"], int(i["START"]), int(i["END"])) for i in isf if i["SEQUENCE"] == c and int(i["END"]) >= st - WIN and int(i["START"]) <= en + WIN]
            mediator = next((m for m, names in MEDIATORS.items() if any(g.startswith(names) for g, _, _ in near)), "none")
            if loc in pt: pt[loc] += 1
            if loc == "chromosome": pv_chrom.add(variant(h["GENE"])); med.add(mediator)
            hit_rows.append({"sample": s, "gene": h["GENE"], "variant": variant(h["GENE"]), "contig": c, "contig_len": L.get(c, 0), "start": st, "end": en,
                             "identity": h["%IDENTITY"], "coverage": h["%COVERAGE"], "location": loc, "platon_plasmid": int(pl is not None and c in pl), "mediator": mediator})
        pc, pp = paper_counts(r["Chromosomal CTX-M variant"]), paper_counts(r["Plasmidic CTX-M variant"])
        oc = sum(h["location"] == "chromosome" for h in ours[s]); op = sum(h["location"] == "plasmid" for h in ours[s])
        om = {h["mediator"] for h in ours[s] if h["location"] == "chromosome"} or {"none"}
        fim = next((re.search(r"FimH(\d+)", i["GENE"], re.I).group(1) for i in abricate(rp / "abricate" / f"{s}.fimtyper.tsv") if re.search(r"FimH(\d+)", i["GENE"], re.I)), "ND")
        ecoh = ",".join(sorted({i["GENE"].split("-")[0] for i in abricate(rp / "abricate" / f"{s}.ecoh.tsv")})) or "ND"
        rows.append({"sample": s, "stratum": r["stratum"],
                     "paper_chrom": sum(pc.values()), "ours_chrom": oc, "tools_chrom": pt["chromosome"],
                     "paper_plasmid": sum(pp.values()), "ours_plasmid": op, "tools_plasmid": pt["plasmid"],
                     "ours_agrees": int(oc == sum(pc.values()) and op == sum(pp.values())),
                     "tools_agrees": int(pt["chromosome"] == sum(pc.values()) and pt["plasmid"] == sum(pp.values())),
                     "paper_chrom_variants": ";".join(sorted(pc)), "tools_chrom_variants": ";".join(sorted(pv_chrom)),
                     "ours_mediator": ",".join(sorted(om)), "tools_mediator": ",".join(sorted(med)) or "none",
                     "platon_large_contig_review": review,
                     "paper_st131_clade": r["Sequence Type 131 clade"], "ours_st131_subclade": typ.get(s, {}).get("st131_subclade", "NA"),
                     "tools_fimH": f"fimH{fim}" if fim != "ND" else "ND", "ours_fimH": typ.get(s, {}).get("fimh", "ND"), "ecoh_otype": ecoh, "ectyper_otype": typ.get(s, {}).get("o_type", "ND"),
                     "paper_serotype": r["Serotype"]})
    for name, data in (("ctxm_hits_paper_tools.tsv", hit_rows), ("per_genome_comparison.tsv", rows)):
        with open(out / name, "w", newline="") as f:
            w = csv.DictWriter(f, list(data[0]), delimiter="\t"); w.writeheader(); w.writerows(data)

    n = len(rows); L = [f"genomes: {n}",
        f"blaCTX-M chromosome+plasmid copy numbers match the paper table: ours {sum(r['ours_agrees'] for r in rows)}/{n}, paper tools {sum(r['tools_agrees'] for r in rows)}/{n}",
        f"chromosomal copies: paper {sum(r['paper_chrom'] for r in rows)}, ours {sum(r['ours_chrom'] for r in rows)}, paper tools {sum(r['tools_chrom'] for r in rows)}",
        f"plasmid copies: paper {sum(r['paper_plasmid'] for r in rows)}, ours {sum(r['ours_plasmid'] for r in rows)}, paper tools {sum(r['tools_plasmid'] for r in rows)}",
        f"ours vs paper tools: same chromosomal count in {sum(r['ours_chrom'] == r['tools_chrom'] for r in rows)}/{n}; same mediator (chromosomal) in {sum(r['ours_mediator'] == r['tools_mediator'] for r in rows)}/{n}",
        f"Platon-called plasmid on a >500 kb contig (flagged for review): {sum(r['platon_large_contig_review'] for r in rows)}"]
    s131 = [r for r in rows if r["paper_st131_clade"] != "non_ST131"]
    L.append(f"ST131 genomes with a paper clade: {len(s131)}; ours matches {sum(r['paper_st131_clade'] == r['ours_st131_subclade'] for r in s131)}")
    L.append(f"fimH ours vs paper-tools (FimTyper via ABRicate): same in {sum(r['ours_fimH'] == r['tools_fimH'] for r in rows)}/{n}")
    mism = [f"{r['sample']}: paper {r['paper_chrom']}c/{r['paper_plasmid']}p ours {r['ours_chrom']}c/{r['ours_plasmid']}p tools {r['tools_chrom']}c/{r['tools_plasmid']}p"
            for r in rows if not (r["ours_agrees"] and r["tools_agrees"])]
    L += ["", "genomes where either pipeline differs from the paper table:"] + mism

    if len(sys.argv) > 4 and Path(sys.argv[4]).exists():
        from Bio import Phylo
        ph = Path(res.parent / "phylo"); trees = {"FastTree (PPanGGOLiN)": sys.argv[4], "IQ-TREE (Panaroo)": ph / "iqtree.treefile", "Mashtree": ph / "mashtree.dnd"}
        T = {}
        for k, p in trees.items():
            if Path(p).exists():
                t = Phylo.read(str(p), "newick")
                for c in t.get_terminals(): c.name = c.name.strip("'")
                T[k] = bipartitions(t)
        with open(out / "tree_rf.tsv", "w") as f:
            f.write("tree_a\ttree_b\tnormalised_RF\n")
            ks = list(T)
            for i, a in enumerate(ks):
                for b in ks[i + 1:]:
                    rf = len(T[a] ^ T[b]) / max(1, len(T[a]) + len(T[b])); f.write(f"{a}\t{b}\t{rf:.3f}\n"); L.append(f"RF {a} vs {b}: {rf:.3f}")
    (out / "summary.txt").write_text("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
