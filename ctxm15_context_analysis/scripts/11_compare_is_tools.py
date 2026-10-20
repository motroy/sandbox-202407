#!/usr/bin/env python3
"""Compare IS-mediator calls for each blaCTX-M hit across four routes:
  mobrecon   MOB-recon mge.report (as used in step 02; IS6-family counted as IS26)
  isfinder   ABRicate against ISfinder (ISEc9 counted as ISEcp1: identical sequence in ISfinder)
  isescan    ISEScan family: IS1380 -> ISEcp1, IS6 -> IS26
  isescan+   ISEScan elements named by BLAST against ISfinder (>= 95% identity over >= 80% of the element)
IS within +-window_bp of the gene count; ISEcp1 takes precedence over IS26 (same as step 02).
Usage: 11_compare_is_tools.py <config.toml> <run_paper_dir> <isescan_windows_dir> <isfinder_db_fasta>   (blastn on PATH)
Writes <run_paper_dir>/comparison/is_mediators.tsv and is_summary.txt.
"""
import csv, subprocess, sys, tomllib
from collections import Counter, defaultdict
from pathlib import Path

ALIAS = {"ISEcp1": "ISEcp1", "ISEc9": "ISEcp1", "IS26": "IS26"}
FAMILY = {"IS1380": "ISEcp1", "IS6": "IS26"}


def pick(names):
    return "ISEcp1" if "ISEcp1" in names else "IS26" if "IS26" in names else "none"


def main():
    cfg = tomllib.loads(Path(sys.argv[1]).read_text()); rp = Path(sys.argv[2]); iw = Path(sys.argv[3]); db = sys.argv[4]
    win = cfg["analysis"].get("window_bp", 10000)
    res = Path(cfg["input"]["outdir"]); out = rp / "comparison"; out.mkdir(exist_ok=True)
    hits = list(csv.DictReader(open(res / "context_hits.tsv"), delimiter="\t"))

    isf, esc, named = defaultdict(list), defaultdict(list), defaultdict(list)
    for s in {h["sample"] for h in hits}:
        p = rp / "abricate" / f"{s}.isfinder.tsv"
        if p.exists():
            for r in csv.DictReader(open(p), delimiter="\t"):
                isf[(s, r["SEQUENCE"])].append((int(r["START"]), int(r["END"]), ALIAS.get(r["GENE"].split("_")[0], "other")))
        tsv = next(iter((iw / s).rglob("*.tsv")), None) if (iw / s).exists() else None
        if tsv is None: continue
        keyed = {}
        for r in csv.DictReader(open(tsv), delimiter="\t"):
            contig, off = r["seqID"].rsplit("@", 1)
            b, e = int(r["isBegin"]) + int(off), int(r["isEnd"]) + int(off)
            esc[(s, contig)].append((min(b, e), max(b, e), FAMILY.get(r["family"], "other")))
            keyed[f"{r['seqID']}_{r['isBegin']}_{r['isEnd']}"] = (contig, min(b, e), max(b, e))
        fna = next(iter(tsv.parent.glob("*.is.fna")), None)
        if fna and keyed:
            q = subprocess.run(["blastn", "-query", str(fna), "-subject", db, "-outfmt", "6 qseqid sseqid pident length qlen", "-max_target_seqs", "1000000"],
                               capture_output=True, text=True).stdout
            best = {}
            for line in q.splitlines():
                qid, sid, pid, ln, ql = line.split("\t"); pid, ln, ql = float(pid), int(ln), int(ql)
                if pid >= 95 and ln >= 0.8 * ql and (qid not in best or (pid, ln) > best[qid][:2]):
                    best[qid] = (pid, ln, ALIAS.get(sid.split("_")[0], "other"))
            for h in (l[1:].split()[0] for l in open(fna) if l.startswith(">")):
                key = h.rsplit("_", 1)[0]
                if key in keyed:
                    c, b, e = keyed[key]; named[(s, c)].append((b, e, best.get(h, (0, 0, "other"))[2]))

    rows = []
    for h in hits:
        s, c, st, en = h["sample"], h["contig"], int(h["start"]), int(h["end"])
        near = lambda d: {n for b, e, n in d.get((s, c), []) if e >= st - win and b <= en + win}
        rows.append({"sample": s, "gene": h["gene"], "location": h["location"], "contig": c, "start": st, "end": en,
                     "mobrecon": h["mediator"], "isfinder": pick(near(isf)), "isescan": pick(near(esc)), "isescan+": pick(near(named))})
    with open(out / "is_mediators.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

    methods = ["mobrecon", "isfinder", "isescan", "isescan+"]
    ch = [r for r in rows if r["location"] == "chromosome"]
    L = [f"chromosomal blaCTX-M copies: {len(ch)}  (paper, all 2,269 genomes: ISEcp1 77.3%, IS26 18.9%)"]
    for m in methods:
        c = Counter(r[m] for r in ch)
        L.append(f"  {m:10s} ISEcp1 {c['ISEcp1']:3d} ({100 * c['ISEcp1'] / len(ch):.0f}%)  IS26 {c['IS26']:3d} ({100 * c['IS26'] / len(ch):.0f}%)  none {c['none']:3d}")
    L.append("pairwise agreement (same call) on chromosomal copies:")
    for i, a in enumerate(methods):
        for b in methods[i + 1:]:
            L.append(f"  {a} vs {b}: {sum(r[a] == r[b] for r in ch)}/{len(ch)}")
    (out / "is_summary.txt").write_text("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
