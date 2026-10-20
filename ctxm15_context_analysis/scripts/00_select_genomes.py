#!/usr/bin/env python3
"""Genome collection: NCBI metadata filter, blaCTX-M detection filter, duplicate removal.

  filter <metadata.tsv> <out_candidates.tsv>
      taxid 562, ANI "OK", N50 >= min_n50, completeness >= min, contamination <= max, BioProject and
      release date present, long-read / hybrid technology. (FASTA + GBFF availability is checked by
      01_download step, or add the files to the table.)
  dedup <candidates.tsv> <resfinder_dir> <cgmlst_dists.tsv> <out_final.tsv>
      keeps genomes with >=1 blaCTX-M (ResFinder identity & coverage >= 80%), then merges genomes into
      one cluster when ALL hold: same BioProject, cgMLST distance <= max_cgmlst_dist, identical ARG
      profile. One representative per cluster: highest N50, then completeness, then accession
      (the tie-break is my assumption; the Methods do not state it).
<resfinder_dir>/<accession>.tsv = ResFinder tab output (columns 'Resistance gene', 'Identity', 'Coverage').
<cgmlst_dists.tsv> = cgmlst-dists output (square matrix, first cell = header label).
Standard library only.
"""
import csv, re, sys, tomllib
from itertools import combinations
from pathlib import Path


def num(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def do_filter(cfg, meta, out):
    c = cfg["selection"]
    keep, drop = [], 0
    rows = list(csv.DictReader(open(meta), delimiter="\t"))
    for r in rows:
        n50, comp, cont = num(r.get(c["n50_col"])), num(r.get(c["completeness_col"])), num(r.get(c["contamination_col"]))
        ok = (r.get(c["taxid_col"], "").strip() == c["taxid"]
              and r.get(c["ani_col"], "").strip() == "OK"
              and n50 is not None and n50 >= c["min_n50"]
              and comp is not None and comp >= c["min_completeness"]
              and cont is not None and cont <= c["max_contamination"]
              and r.get(c["bioproject_col"], "").strip()
              and r.get(c["release_date_col"], "").strip()
              and re.search(c["long_read_regex"], r.get(c["tech_col"], "")))
        keep.append(r) if ok else None
        drop += not ok
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, list(rows[0]) if rows else [], delimiter="\t"); w.writeheader(); w.writerows(keep)
    print(f"{len(keep)} kept, {drop} excluded of {len(rows)}", file=sys.stderr)


def resfinder_profile(path, min_id, min_cov):
    genes = set()
    for r in csv.DictReader(open(path), delimiter="\t"):
        i, cv = num(r.get("Identity")), num(r.get("Coverage"))
        if i is not None and cv is not None and i >= min_id and cv >= min_cov:
            genes.add(r["Resistance gene"])
    return frozenset(genes)


def do_dedup(cfg, cand, rdir, dist, out):
    c = cfg["selection"]
    acc = c["accession_col"]
    rows = {r[acc]: r for r in csv.DictReader(open(cand), delimiter="\t")}
    prof = {}
    for a in rows:
        p = Path(rdir) / f"{a}.tsv"
        if p.exists():
            g = resfinder_profile(p, c["resfinder_min_identity"], c["resfinder_min_coverage"])
            if any(x.startswith("blaCTX-M") for x in g):
                prof[a] = g
    print(f"{len(prof)} of {len(rows)} candidates carry blaCTX-M", file=sys.stderr)

    with open(dist) as f:
        rd = csv.reader(f, delimiter="\t"); cols = next(rd)[1:]
        D = {r[0]: dict(zip(cols, map(float, r[1:]))) for r in rd}
    parent = {a: a for a in prof}
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for a, b in combinations(prof, 2):
        if (rows[a][c["bioproject_col"]] == rows[b][c["bioproject_col"]] and prof[a] == prof[b]
                and a in D and b in D[a] and D[a][b] <= c["max_cgmlst_dist"]):
            parent[find(a)] = find(b)
    clusters = {}
    for a in prof: clusters.setdefault(find(a), []).append(a)
    key = lambda a: (-(num(rows[a][c["n50_col"]]) or 0), -(num(rows[a][c["completeness_col"]]) or 0), a)
    final = [min(m, key=key) for m in clusters.values()]
    with open(out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(list(next(iter(rows.values()))) + ["cluster_size", "arg_profile"])
        for a in sorted(final):
            w.writerow([*rows[a].values(), len(clusters[find(a)]), ",".join(sorted(prof[a]))])
    print(f"{len(final)} non-redundant genomes", file=sys.stderr)


if __name__ == "__main__":
    cfg = tomllib.loads(Path("config.toml").read_text())
    if len(sys.argv) == 4 and sys.argv[1] == "filter": do_filter(cfg, *sys.argv[2:])
    elif len(sys.argv) == 6 and sys.argv[1] == "dedup": do_dedup(cfg, *sys.argv[2:])
    else: sys.exit(__doc__)
