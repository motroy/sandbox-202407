"""FWSEC0007: assembly accuracy vs truth + AMR call matrix (truth / assembly stages / Amira). usage: fwsec_summary.py OUTDIR"""
import math, re, subprocess, sys
from pathlib import Path
import pandas as pd

R = Path("results/FWSEC0007"); out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
truth = "data/FWSEC0007/truth.fna"
stages = {"myloasm": R/"polish/myloasm.fasta", "polypolish": R/"polish/polypolish.fasta", "polished": R/"polish/polypolish_pypolca.fasta"}

# --- assembly accuracy: minimap2 asm5 alignments to the complete genome (primary alignments)
acc = []
for name, fa in stages.items():
    paf = subprocess.run(["minimap2", "-cx", "asm5", "--secondary=no", truth, str(fa)], capture_output=True, text=True, check=True).stdout
    nm = blk = 0
    for l in paf.splitlines():
        f = l.split("\t")
        tags = {t.split(":")[0]: t.split(":")[2] for t in f[12:]}
        nm += int(tags["NM"]); blk += int(f[10])
    lens = [len(s) for s in "".join(open(fa)).split(">")[1:] for s in ["".join(s.split("\n")[1:])]]
    acc.append({"stage": name, "contigs": len(lens), "total_bp": sum(lens), "aligned_bp": blk, "edits(NM)": nm,
                "error_rate_%": round(100*nm/blk, 4), "QV": round(-10*math.log10(max(nm,1)/blk), 1)})
pd.DataFrame(acc).to_csv(out/"assembly_accuracy.tsv", sep="\t", index=False)

# --- AMR call matrix
def fam(n): return re.sub(r"^bla(CTX-M|TEM|OXA|SHV|KPC|CMY|NDM|EC)-?.*$", r"bla\1", n)
cols = {}
for k in ["truth", "myloasm", "polypolish", "polished"]:
    d = pd.read_csv(R/f"amrfinder/{k}.tsv", sep="\t"); d = d[d["Type"] == "AMR"]
    cols[k] = d.apply(lambda r: ("[" + r["Method"] + "] " if r["Subtype"] != "AMR" or r["Method"] not in ("BLASTX","EXACTX","ALLELEX","COMPLETE") else "") + fam(r["Element symbol"]), axis=1).tolist()
a = pd.read_csv(R/"amira/amira_results.tsv", sep="\t")
cols["amira"] = sorted({fam(x.split("/")[0]) for x in a["Determinant name"]})
allg = sorted({g for v in cols.values() for g in v})
m = pd.DataFrame({k: [v.count(g) for g in allg] for k, v in cols.items()}, index=allg)
m.to_csv(out/"amr_call_matrix.tsv", sep="\t")
print(pd.DataFrame(acc).to_string(index=False)); print(); print(m.to_string())
