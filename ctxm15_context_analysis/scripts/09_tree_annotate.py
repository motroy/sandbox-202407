#!/usr/bin/env python3
"""Annotate the Mashtree / IQ-TREE phylogenies with typing and blaCTX-M integration status.

Usage: 09_tree_annotate.py <config.toml> <phylo_dir>      (phylo_dir from 08_phylogeny.sh; needs results/ from 02 and 05)
Writes <phylo_dir>/annotation/:
  tree_metadata.tsv      per genome: ST, phylogroup, ST131 subclade, blaCTX-M variants, integration status, mediators
  itol_<field>.txt       iTOL colour-strip datasets (st, integration, variant)  - upload with the .treefile / .dnd
  tree_<name>.png        rectangular tree (midpoint-rooted) with coloured strips
  tree_comparison.txt    normalised Robinson-Foulds distance between the Mashtree and IQ-TREE topologies
Needs Biopython and matplotlib.
"""
import csv, sys, tomllib
from pathlib import Path

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

PAL = ["#8dd3c7", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#a6cee3"]
INTEG = {"chromosome": "#f28e8e", "both": "#c9a9e3", "plasmid": "#8fbbe8", "none": "#eeeeee"}


def load_meta(res, samples):
    typ = {r["sample"]: r for r in csv.DictReader(open(res / "typing.tsv"), delimiter="\t")} if (res / "typing.tsv").exists() else {}
    hits = {}
    for h in csv.DictReader(open(res / "context_hits.tsv"), delimiter="\t"):
        hits.setdefault(h["sample"], []).append(h)
    meta = {}
    for s in samples:
        hs = hits.get(s, []); locs = {h["location"] for h in hs}
        t = typ.get(s, {})
        meta[s] = {"sample": s, "st": t.get("st", "ND"), "phylogroup": t.get("phylogroup", "ND"), "st131_subclade": t.get("st131_subclade", "NA"),
                   "variant": ";".join(sorted({h["gene"] for h in hs})) or "none",
                   "integration": "both" if {"chromosome", "plasmid"} <= locs else "chromosome" if "chromosome" in locs else "plasmid" if "plasmid" in locs else "none",
                   "mediators": ",".join(sorted({h["mediator"] for h in hs if h["location"] == "chromosome"})) or "none"}
    return meta


def match(tip, samples):
    """Tree tip -> sample id (Panaroo drops the final '.1' of file names)."""
    if tip in samples: return tip
    c = [s for s in samples if s.startswith(tip + ".") or s.replace(".", "_") == tip]
    return c[0] if len(c) == 1 else None


def bipartitions(tree):
    tips = {t.name for t in tree.get_terminals()}
    out = set()
    for cl in tree.get_nonterminals():
        s = frozenset(t.name for t in cl.get_terminals())
        if 1 < len(s) < len(tips) - 1: out.add(min(s, tips - s, key=lambda x: sorted(x)))
    return out


def colours(values):
    uniq = sorted(set(values), key=lambda v: (-list(values).count(v), v))
    return {v: PAL[i % len(PAL)] for i, v in enumerate(uniq)}


def itol(path, field, cmap, ids):
    with open(path, "w") as f:
        f.write(f"DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\t{field}\nCOLOR\t#000000\nLEGEND_TITLE\t{field}\n")
        f.write("LEGEND_SHAPES\t" + "\t".join("1" for _ in cmap) + "\nLEGEND_COLORS\t" + "\t".join(cmap.values()) + "\nLEGEND_LABELS\t" + "\t".join(cmap) + "\nDATA\n")
        for tip, (val) in ids.items(): f.write(f"{tip}\t{cmap[val]}\t{val}\n")


def draw(tree, name, meta, samples, out):
    tree.root_at_midpoint(); tree.ladderize()
    tips = tree.get_terminals(); depth = tree.depths()
    ypos = {t: i for i, t in enumerate(tips)}
    for cl in tree.get_nonterminals(order="postorder"): ypos[cl] = sum(ypos[c] for c in cl.clades) / len(cl.clades)
    fig, ax = plt.subplots(figsize=(11, 0.28 * len(tips) + 2))
    for cl in tree.find_clades(order="level"):
        x0 = depth[cl]
        for c in cl.clades:
            ax.plot([x0, depth[c]], [ypos[c], ypos[c]], "k-", lw=0.9)
        if cl.clades: ax.plot([x0, x0], [ypos[cl.clades[0]], ypos[cl.clades[-1]]], "k-", lw=0.9)
    xmax = max(depth.values()); w = xmax * 0.05
    fields = [("integration", INTEG), ("st", colours([m["st"] for m in meta.values()])),
              ("variant", colours([m["variant"] for m in meta.values()]))]
    for k, (field, cmap) in enumerate(fields):
        for t in tips:
            s = match(t.name, samples)
            if s: ax.add_patch(plt.Rectangle((xmax * 1.03 + k * w * 1.15, ypos[t] - 0.42), w, 0.84, color=cmap[meta[s][field]], lw=0))
        ax.text(xmax * 1.03 + k * w * 1.15 + w / 2, -1.4, field, rotation=60, ha="left", va="bottom", fontsize=8)
    for t in tips:
        s = match(t.name, samples)
        ax.text(xmax * 1.03 + 3 * w * 1.15 + w * 0.3, ypos[t], f"{t.name}  {meta[s]['st'] if s else ''}", va="center", fontsize=6.5)
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in INTEG.values()]
    ax.legend(handles, [f"{k} blaCTX-M" for k in INTEG], loc="lower left", fontsize=7, title="integration", title_fontsize=7, frameon=False)
    ax.set_ylim(len(tips), -2.5); ax.set_xlim(-xmax * 0.02, xmax * 2.1); ax.axis("off"); ax.set_title(name, fontsize=10, loc="left")
    fig.savefig(out, dpi=170, bbox_inches="tight", facecolor="white"); plt.close(fig)


def main():
    cfg = tomllib.loads(Path(sys.argv[1]).read_text()); phylo = Path(sys.argv[2])
    res = Path(cfg["input"]["outdir"]); ann = phylo / "annotation"; ann.mkdir(exist_ok=True)
    samples = [r["sample"] for r in csv.DictReader(open(cfg["input"]["samples"]), delimiter="\t")]
    meta = load_meta(res, samples)
    with open(ann / "tree_metadata.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, list(next(iter(meta.values()))), delimiter="\t"); w.writeheader(); w.writerows(meta.values())
    trees = {"IQ-TREE (Panaroo core alignment)": phylo / "iqtree.treefile", "Mashtree": phylo / "mashtree.dnd"}
    read = {}
    for name, p in trees.items():
        if not p.exists(): continue
        t = Phylo.read(str(p), "newick")
        for c in t.get_terminals(): c.name = (c.name or "").strip("'")
        read[name] = t
        stem = "iqtree" if name.startswith("IQ") else "mashtree"
        ids = {c.name: match(c.name, samples) for c in t.get_terminals()}
        for field in ("st", "integration", "variant"):
            cm = INTEG if field == "integration" else colours([m[field] for m in meta.values()])
            itol(ann / f"itol_{field}_{stem}.txt", field, cm, {tip: meta[s][field] for tip, s in ids.items() if s})
        draw(t, name, meta, samples, ann / f"tree_{stem}.png")
    if len(read) == 2:
        a, b = list(read.values()); shared = {t.name for t in a.get_terminals()} & {t.name for t in b.get_terminals()}
        for t in (a, b): t.prune  # topology compared on shared tips only
        if len(shared) < len(a.get_terminals()):
            for t in (a, b):
                for c in [c for c in t.get_terminals() if c.name not in shared]: t.prune(c)
        ba, bb = bipartitions(a), bipartitions(b); rf = len(ba ^ bb) / max(1, len(ba) + len(bb))
        (ann / "tree_comparison.txt").write_text(f"tips compared: {len(shared)}\nnormalised Robinson-Foulds (0 = identical, 1 = no shared splits): {rf:.3f}\n")
        print(open(ann / "tree_comparison.txt").read(), file=sys.stderr)


if __name__ == "__main__":
    main()
