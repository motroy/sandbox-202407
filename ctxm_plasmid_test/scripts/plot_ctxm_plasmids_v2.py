#!/usr/bin/env python3
"""
Reference-guided backbone alignment of blaCTX-M (or any --gene) plasmids from SHORT-READ (fragmented) assemblies.

Per MOB-suite cluster (e.g. AA147):
  1. Each isolate's plasmid_<CLUSTER>.fasta (ALL contigs MOB-recon assigned to that cluster) is
     treated as one plasmid "bin".
  2. A reference is chosen: the bin with the longest single contig in the cluster (any isolate,
     CTX-M or not), or one you supply with --refs AA147=ref.fasta.
  3. Every bin's contigs are blasted against the reference, then ordered and oriented by their
     position on it, and joined with a visible gap (white box) = contig break.
     Contigs with no match to the reference are appended at the end ("unplaced").
  4. Tracks are sorted by how much of the reference they cover; adjacent tracks are linked by
     blastn ribbons (grey = same strand, red = inverted).
  5. AMRFinderPlus genes are drawn as arrows (colour = drug class), the target gene (--gene, default
     blaCTX-M = all variants; use blaCTX-M-15 for that allele only) outlined in black.
     A track label ends with "*" when the blaCTX-M contig itself could NOT be placed on the
     reference (so it cannot be claimed to sit on that backbone).

Expected input:
  <mobrecon_dir>/<sample>/plasmid_<CLUSTER>.fasta
  <amr_dir>/<sample>*.tsv   (AMRFinderPlus; "Contig id" == first word of the FASTA header)
  optional --metadata TSV with columns: sample, host

Requires: pandas biopython pygenomeviz matplotlib ; blastn in PATH
"""
import argparse
import tempfile
import subprocess
try:
    import tomllib            # Python >= 3.11
except ModuleNotFoundError:   # pragma: no cover
    import tomli as tomllib
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pygenomeviz import GenomeViz

CTX = "blaCTX-M"   # default target gene; override with --gene (e.g. blaCTX-M-15)
GAP = 500  # bp drawn between contigs
DEFAULT_CLUSTERS = ["AA147", "AC315", "AA474", "AA161", "AA176", "AB690", "AA887", "AA738"]


def is_target(symbol, gene):
    """True for the target gene or any of its variants: gene 'blaCTX-M' matches blaCTX-M and blaCTX-M-15 etc.,
    gene 'blaCTX-M-15' matches only blaCTX-M-15 (not blaCTX-M-1 or blaCTX-M-150). Case-insensitive."""
    s, g = str(symbol).lower(), gene.lower()
    return s == g or s.startswith(g + "-")


AMR_COLS = ["contig", "start", "end", "strand", "symbol", "cls"]


# ------------------------------------------------------------------ input
def read_amr(path):
    df = pd.read_csv(path, sep="\t")
    ren = {"Gene symbol": "symbol", "Element symbol": "symbol", "Contig id": "contig",
           "Start": "start", "Stop": "end", "Strand": "strand", "Class": "cls"}
    df = df.rename(columns={k: v for k, v in ren.items() if k in df.columns})
    if "cls" not in df.columns:
        df["cls"] = "OTHER"
    df["cls"] = df["cls"].fillna("OTHER")
    return df[AMR_COLS].copy()


def find_amr_file(amr_dir, sample):
    for pat in (f"{sample}.amrfinder.tsv", f"{sample}*.tsv", f"{sample}*amr*.tsv", f"{sample}/*.tsv"):
        hits = sorted(amr_dir.glob(pat))
        if hits:
            return hits[0]
    return None


MGE_COLS = ["contig", "start", "end", "strand", "name"]


def read_mge(path):
    """Parse MOB-recon mge.report.txt. Column names are matched loosely; if your header differs
    the error message prints it so the matching can be adjusted."""
    if path is None or not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=MGE_COLS)
    df = pd.read_csv(path, sep="\t")
    low = {c: c.lower() for c in df.columns}

    def pick(*keys, exclude=()):
        for k in keys:   # keys are in priority order
            for c, l in low.items():
                if k in l and not any(x in l for x in exclude):
                    return c
        return None

    # MOB-suite 3.x mge.report.txt has mge_start/mge_end (position within the MGE reference) AND
    # contig_start/contig_end (position on the contig): use the contig ones, and mge_type (e.g. ISSen6) as label.
    cols = dict(contig=pick("contig_id", "contig"),
                start=pick("contig_start") or pick("start", exclude=("mge_",)) or pick("start"),
                end=pick("contig_end") or pick("end", "stop", exclude=("mge_",)) or pick("end", "stop"),
                strand=pick("strand"),
                name=pick("mge_type", "mge_acs", "element", "family", "type", "name", "subject",
                          exclude=("molecule", "contig", "sample")) or pick("mge_id", "mge"))
    if not all(cols[k] for k in ("contig", "start", "end")):
        raise SystemExit(f"Cannot parse {path}: columns are {list(df.columns)}")
    out = pd.DataFrame({
        "contig": df[cols["contig"]].astype(str).str.split().str[0],
        "start": df[cols["start"]], "end": df[cols["end"]],
        "strand": df[cols["strand"]] if cols["strand"] else "+",
        "name": df[cols["name"]].astype(str) if cols["name"] else "MGE"})
    lo, hi = out[["start", "end"]].min(axis=1), out[["start", "end"]].max(axis=1)
    out["start"], out["end"] = lo.astype(int), hi.astype(int)
    out["strand"] = out["strand"].map(lambda x: "-" if str(x) in ("-", "-1", "minus") else "+")
    return out[MGE_COLS]


def load_bins(mob_dir, amr_dir, clusters, meta, mge_name="mge.report.txt", gene=CTX):
    bins = defaultdict(list)
    for sd in sorted(p for p in mob_dir.iterdir() if p.is_dir()):
        amr_file = find_amr_file(amr_dir, sd.name)
        amr = read_amr(amr_file) if amr_file else pd.DataFrame(columns=AMR_COLS)
        mge_all = read_mge(sd / mge_name)
        for cl in clusters:
            fa = sd / f"plasmid_{cl}.fasta"
            if not fa.exists():
                continue
            contigs = list(SeqIO.parse(fa, "fasta"))
            ids = {c.id for c in contigs}
            hits = amr[amr["contig"].isin(ids)].reset_index(drop=True)
            ctx = hits["symbol"].map(lambda s: is_target(s, gene))
            bins[cl].append(dict(sample=sd.name, host=meta.get(sd.name, ""), contigs=contigs, amr=hits,
                                 mge=mge_all[mge_all["contig"].isin(ids)].reset_index(drop=True),
                                 has_ctx=bool(ctx.any()),
                                 ctx_contigs=set(hits.loc[ctx, "contig"])))
    return bins


# ------------------------------------------------------------------ blast + assembly
def blast(q_recs, s_recs, min_id, min_len):
    with tempfile.TemporaryDirectory() as td:
        q, s = Path(td, "q.fa"), Path(td, "s.fa")
        SeqIO.write(q_recs, q, "fasta")
        SeqIO.write(s_recs, s, "fasta")
        out = subprocess.run(
            ["blastn", "-query", str(q), "-subject", str(s), "-task", "megablast", "-evalue", "1e-10",
             "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send"],
            capture_output=True, text=True, check=True).stdout
    rows = []
    for line in out.strip().splitlines():
        qid, sid, pid, aln, qs, qe, ss, se = line.split("\t")
        if float(pid) >= min_id and int(aln) >= min_len:
            rows.append((qid, sid, float(pid), int(aln), int(qs), int(qe), int(ss), int(se)))
    return rows


def assemble(b, ref_rec, min_id, min_len):
    """Order/orient the contigs of bin b along ref_rec and concatenate them."""
    info = {}
    if ref_rec is not None:
        for qid, _, pid, aln, qs, qe, ss, se in blast(b["contigs"], [ref_rec], min_id, min_len):
            d = info.setdefault(qid, dict(plus=0, minus=0, pos=0.0, w=0))
            d["plus" if se >= ss else "minus"] += aln
            d["pos"] += min(ss, se) * aln
            d["w"] += aln
    placed = sorted(((d["pos"] / d["w"], d["minus"] > d["plus"], c)
                     for c in b["contigs"] if (d := info.get(c.id))), key=lambda x: x[0])
    unplaced = sorted((c for c in b["contigs"] if c.id not in info), key=lambda c: -len(c.seq))
    order = [(flip, c) for _, flip, c in placed] + [(False, c) for c in unplaced]

    pieces, genes, mges, breaks, off = [], [], [], [], 0
    for k, (flip, c) in enumerate(order):
        L = len(c.seq)
        seq = c.seq.reverse_complement() if flip else c.seq
        for _, m in b["mge"][b["mge"]["contig"] == c.id].iterrows():
            st, en, sd = max(int(m["start"]), 1), min(int(m["end"]), L), m["strand"]   # MOB-recon can overhang contig ends
            if en < st:
                continue
            if flip:
                st, en, sd = L - en + 1, L - st + 1, ("-" if sd == "+" else "+")
            mges.append(dict(name=m["name"], start=st + off, end=en + off, strand=sd))
        for _, g in b["amr"][b["amr"]["contig"] == c.id].iterrows():
            st, en, sd = max(int(g["start"]), 1), min(int(g["end"]), L), g["strand"]
            if en < st:
                continue
            if flip:
                st, en, sd = L - en + 1, L - st + 1, ("-" if sd == "+" else "+")
            genes.append(dict(symbol=g["symbol"], start=st + off, end=en + off, strand=sd, cls=g["cls"]))
        pieces.append(str(seq))
        off += L
        if k < len(order) - 1:
            pieces.append("N" * GAP)
            breaks.append((off, off + GAP))
            off += GAP
    rec = SeqRecord(Seq("".join(pieces)), id=b["sample"], description="")
    b2 = dict(b)
    b2.update(rec=rec, genes=pd.DataFrame(genes, columns=["symbol", "start", "end", "strand", "cls"]),
              mges=pd.DataFrame(mges, columns=["name", "start", "end", "strand"]),
              breaks=breaks, n_unplaced=len(unplaced), total=sum(len(c.seq) for c in b["contigs"]),
              aligned_bp=sum(d["w"] for d in info.values()),
              ctx_placed=all(cid in info for cid in b["ctx_contigs"]) if b["ctx_contigs"] else None)
    return b2


# ------------------------------------------------------------------ plotting
def stagger(items, maxlen, min_gap=0.024):
    """Label placement without overlap. items: list of (start_bp, label). Returns a list of
    (vpos, raise_lines) per item. Labels are rotated 45 degrees, so two labels on the same level only collide
    when their start points are closer than about one text line (min_gap * longest track, in bp). Levels:
    above the track, above and raised by 4 and 8 blank lines, below the track. Greedy left to right."""
    levels = [("top", 0), ("top", 4), ("top", 8), ("bottom", 0)]
    last = [-1e18] * len(levels)
    out = [None] * len(items)
    for i in sorted(range(len(items)), key=lambda k: items[k][0]):
        start = items[i][0]
        lvl = next((k for k in range(len(levels)) if start - last[k] >= min_gap * maxlen), None)
        if lvl is None:                      # all levels busy: use the one that is free the longest
            lvl = min(range(len(levels)), key=lambda k: last[k])
        last[lvl] = start
        out[i] = levels[lvl]
    return out


def draw(cl, tracks, ref_label, outdir, min_id, gene=CTX):
    classes = sorted({c for t in tracks for c in t["genes"]["cls"]})
    cmap = plt.get_cmap("tab20")
    pal = {c: cmap(i % 20) for i, c in enumerate(classes)}

    gv = GenomeViz(track_align_type="left", fig_track_height=0.6)
    maxlen = max(len(t["rec"].seq) for t in tracks)
    names = []
    for i, t in enumerate(tracks):
        flag = "*" if t["ctx_placed"] is False else ""
        name = f'{i}_{t["sample"]}{flag}'
        label = name.split("_", 1)[1] + (f' | {t["host"]}' if t["host"] else "") + \
            f' | {t["total"]/1000:.0f} kb, {len(t["contigs"])} contigs'
        tr = gv.add_feature_track(label, len(t["rec"].seq))
        names.append(label)
        for s, e in t["breaks"]:
            tr.add_feature(s, e, 1, plotstyle="box", fc="white", ec="none")
        # label placement: genes and MGEs share the staggering so neighbouring labels do not overlap
        items = [(int(m["start"]), str(m["name"])) for _, m in t["mges"].iterrows()] + \
                [(int(g["start"]), str(g["symbol"])) for _, g in t["genes"].iterrows()]
        place = stagger(items, maxlen)
        nm = len(t["mges"])

        def txt(k, color):
            vpos, up = place[k]
            return items[k][1] + "\n" * up, dict(rotation=45, vpos=vpos, hpos="left", size=7, color=color)

        for k, (_, m) in enumerate(t["mges"].iterrows()):   # MOB-recon MGEs: dark grey boxes
            lab, kw = txt(k, "#333333")
            tr.add_feature(int(m["start"]), int(m["end"]), 1 if m["strand"] == "+" else -1,
                           plotstyle="box", fc="dimgray", ec="none", label=lab, text_kws=kw)
        for k, (_, g) in enumerate(t["genes"].iterrows(), start=nm):
            is_ctx = is_target(g["symbol"], gene)
            lab, kw = txt(k, "black")
            tr.add_feature(int(g["start"]), int(g["end"]), 1 if g["strand"] == "+" else -1,
                           plotstyle="bigarrow", fc=pal[g["cls"]], ec="black" if is_ctx else "none",
                           lw=1.5 if is_ctx else 0, label=lab, text_kws=kw)

    for (n1, t1), (n2, t2) in zip(zip(names[:-1], tracks[:-1]), zip(names[1:], tracks[1:])):
        a = SeqRecord(t1["rec"].seq, id="a", description="")
        b = SeqRecord(t2["rec"].seq, id="b", description="")
        for _, _, pid, _, qs, qe, ss, se in blast([a], [b], min_id, 300):
            gv.add_link((n1, qs, qe), (n2, ss, se), color="lightgrey", inverted_color="tomato",
                        v=pid, vmin=min_id, curve=False)

    fig = gv.plotfig()
    try:
        gv.set_colorbar(fig, bar_colors=["lightgrey"], vmin=min_id, bar_label="blastn identity (%)")
    except Exception:
        pass
    fig.legend(handles=[Patch(fc=c, label=k) for k, c in pal.items()] + [Patch(fc="dimgray", label="MGE (MOB-recon)")],
               title="AMR class",
               loc="upper center", ncol=min(len(pal), 5), fontsize=8, frameon=False,
               bbox_to_anchor=(0.5, 0.0))
    fig.canvas.draw()   # gene labels stick out above the figure box: put the title above everything drawn so far
    bb = fig.get_tightbbox(fig.canvas.get_renderer())
    top = bb.y1 / fig.get_size_inches()[1] + 0.02
    fig.suptitle(f"{cl}: {gene} plasmids, ordered on reference {ref_label} "
                 f"(white = contig break; * = {gene} contig not placed on reference)", fontsize=9,
                 y=top, va="bottom")
    for ext in ("png", "svg", "pdf"):
        fig.savefig(outdir / f"{cl}_ctxm_alignment.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", type=Path, help="TOML config file (CLI options override it)")
    ap.add_argument("--mobrecon-dir", type=Path)
    ap.add_argument("--amr-dir", type=Path)
    ap.add_argument("--metadata", type=Path)
    ap.add_argument("--gene", default=CTX, help="target gene or variant, e.g. blaCTX-M (all variants) or blaCTX-M-15 (only that allele)")
    ap.add_argument("--clusters", nargs="+", default=DEFAULT_CLUSTERS)
    ap.add_argument("--refs", nargs="*", default=[], help="CLUSTER=ref.fasta (optional external reference)")
    ap.add_argument("--outdir", type=Path, default=Path("ctxm_plasmid_figs"))
    ap.add_argument("--mge-file", default="mge.report.txt", help="MOB-recon MGE report name inside each sample dir")
    ap.add_argument("--max-per-cluster", type=int, default=8)
    ap.add_argument("--min-aln-len", type=int, default=300)
    ap.add_argument("--min-identity", type=float, default=90.0)
    pre, _ = ap.parse_known_args()
    if pre.config:
        with open(pre.config, "rb") as fh:
            cfg = {k.replace("-", "_"): v for k, v in tomllib.load(fh).items()}
        unknown = set(cfg) - {x.dest for x in ap._actions}
        if unknown:
            ap.error(f"unknown key(s) in {pre.config}: {sorted(unknown)}")
        for k in ("mobrecon_dir", "amr_dir", "metadata", "outdir"):
            if k in cfg:
                cfg[k] = Path(cfg[k])
        if isinstance(cfg.get("refs"), dict):          # [refs] table -> CLUSTER=path list
            cfg["refs"] = [f"{k}={v}" for k, v in cfg["refs"].items()]
        ap.set_defaults(**cfg)
    a = ap.parse_args()
    if not (a.mobrecon_dir and a.amr_dir):
        ap.error("mobrecon_dir and amr_dir are required (on the command line or in the config file)")
    a.outdir.mkdir(parents=True, exist_ok=True)

    meta = {}
    if a.metadata:
        m = pd.read_csv(a.metadata, sep="\t")
        meta = dict(zip(m["sample"].astype(str), m["host"].astype(str)))
    ext_refs = dict(r.split("=", 1) for r in a.refs)

    bins = load_bins(a.mobrecon_dir, a.amr_dir, a.clusters, meta, a.mge_file, a.gene)
    summary = []
    for cl in a.clusters:
        allb = bins.get(cl, [])
        ctxb = [b for b in allb if b["has_ctx"]]
        if len(ctxb) < 2:
            print(f"[skip] {cl}: {len(ctxb)} {a.gene} isolate(s)")
            continue

        # reference
        if cl in ext_refs:
            rc = list(SeqIO.parse(ext_refs[cl], "fasta"))
            rb = dict(sample=f"ref_{cl}", host="", contigs=rc, amr=pd.DataFrame(columns=AMR_COLS),
                      mge=pd.DataFrame(columns=MGE_COLS),
                      has_ctx=False, ctx_contigs=set())
        else:
            rb = max(allb, key=lambda b: (max(len(c.seq) for c in b["contigs"]), -len(b["contigs"])))
        ref = assemble(rb, None, a.min_identity, a.min_aln_len)   # reference: contigs by length
        ref["ctx_placed"] = True if ref["ctx_contigs"] else None  # reference defines the coordinates
        ref_rec = SeqRecord(ref["rec"].seq, id="ref", description="")

        asm = [assemble(b, ref_rec, a.min_identity, a.min_aln_len) for b in ctxb if b["sample"] != rb["sample"]]
        asm.sort(key=lambda t: -t["aligned_bp"])
        # cap: one per host first, then by coverage
        picked, seen = [], set()
        for t in asm:
            if t["host"] not in seen:
                picked.append(t); seen.add(t["host"])
        picked += [t for t in asm if t not in picked]
        picked = sorted(picked[:a.max_per_cluster - 1], key=lambda t: -t["aligned_bp"])
        tracks = [ref] + picked

        for t in [ref] + asm:
            ctx = t["genes"][t["genes"]["symbol"].map(lambda s: is_target(s, a.gene))]
            is_ref = t is ref
            summary.append(dict(cluster=cl, sample=t["sample"], host=t["host"], n_contigs=len(t["contigs"]),
                                total_bp=t["total"], aligned_bp_to_ref=t["total"] if is_ref else t["aligned_bp"],
                                n_unplaced_contigs=0 if is_ref else t["n_unplaced"], ctxm=";".join(ctx["symbol"]),
                                ctxm_contig_placed=t["ctx_placed"], reference=rb["sample"], is_reference=is_ref))
        print(f"[plot] {cl}: ref={rb['sample']}, {len(tracks)} tracks")
        draw(cl, tracks, rb["sample"], a.outdir, a.min_identity, a.gene)

    pd.DataFrame(summary).to_csv(a.outdir / "ctxm_plasmid_summary.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
