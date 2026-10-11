"""FWSEC0007: assembly accuracy vs truth + AMR call matrix (truth / assemblers x polishing stage / Amira). usage: fwsec_summary.py OUTDIR"""
import math, re, subprocess, sys
from pathlib import Path
import pandas as pd

R = Path("results/FWSEC0007"); out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
truth = "data/FWSEC0007/truth.fna"
stages = {  # label -> (fasta, amrfinder tsv)
    "myloasm": (R/"polish/myloasm.fasta", R/"amrfinder/myloasm.tsv"),
    "myloasm+polypolish": (R/"polish/polypolish.fasta", R/"amrfinder/polypolish.tsv"),
    "myloasm+polypolish+pypolca": (R/"polish/polypolish_pypolca.fasta", R/"amrfinder/polished.tsv"),
    "dragonflye(flye)": (R/"dragonflye/asm/flye-unpolished.fasta", R/"dragonflye/amrfinder/flye.tsv"),
    "dragonflye(flye+racon)": (R/"dragonflye/polish/input.fasta", R/"dragonflye/amrfinder/racon.tsv"),
    "dragonflye+polypolish": (R/"dragonflye/polish/polypolish.fasta", R/"dragonflye/amrfinder/polypolish.tsv"),
    "dragonflye+polypolish+pypolca": (R/"dragonflye/polish/polypolish_pypolca.fasta", R/"dragonflye/amrfinder/polished.tsv"),
}

def seq_lengths(fa):
    L, cur = [], 0
    for l in open(fa):
        if l.startswith(">"):
            if cur: L.append(cur)
            cur = 0
        else: cur += len(l.strip())
    if cur: L.append(cur)
    return L

tlen = sum(seq_lengths(truth))
acc = []
for name, (fa, _) in stages.items():
    paf = subprocess.run(["minimap2", "-cx", "asm5", "--secondary=no", truth, str(fa)], capture_output=True, text=True, check=True).stdout
    nm = blk = 0; iv = {}
    for l in paf.splitlines():
        f = l.split("\t"); tags = {t.split(":")[0]: t.split(":")[2] for t in f[12:]}
        nm += int(tags["NM"]); blk += int(f[10]); iv.setdefault(f[5], []).append((int(f[7]), int(f[8])))
    cov = 0
    for ivs in iv.values():  # union of target intervals = fraction of truth covered
        ivs.sort(); s, e = ivs[0]
        for a, b in ivs[1:]:
            if a > e: cov += e - s; s, e = a, b
            else: e = max(e, b)
        cov += e - s
    lens = seq_lengths(fa)
    acc.append({"stage": name, "contigs": len(lens), "total_bp": sum(lens), "truth_covered_%": round(100*cov/tlen, 1),
                "aligned_bp": blk, "edits(NM)": nm, "error_rate_%": round(100*nm/blk, 4), "QV": round(-10*math.log10(max(nm, 1)/blk), 1)})
pd.DataFrame(acc).to_csv(out/"assembly_accuracy.tsv", sep="\t", index=False)

def fam(n): return re.sub(r"^bla(CTX-M|TEM|OXA|SHV|KPC|CMY|NDM|EC)-?.*$", r"bla\1", n)
def calls(tsv):
    d = pd.read_csv(tsv, sep="\t"); d = d[d["Type"] == "AMR"]
    return d.apply(lambda r: ("[" + r["Method"] + "] " if r["Subtype"] != "AMR" or r["Method"] not in ("BLASTX", "EXACTX", "ALLELEX", "COMPLETE") else "") + fam(r["Element symbol"]), axis=1).tolist()
cols = {"truth": calls(R/"amrfinder/truth.tsv")}
for k, (_, t) in stages.items(): cols[k] = calls(t)
a = pd.read_csv(R/"amira/amira_results.tsv", sep="\t")
cols["amira"] = sorted({fam(x.split("/")[0]) for x in a["Determinant name"]})
allg = sorted({g for v in cols.values() for g in v})
m = pd.DataFrame({k: [v.count(g) for g in allg] for k, v in cols.items()}, index=allg)
m.to_csv(out/"amr_call_matrix.tsv", sep="\t")
pd.set_option("display.width", 250)
print(pd.DataFrame(acc).to_string(index=False)); print(); print(m.T.to_string())
