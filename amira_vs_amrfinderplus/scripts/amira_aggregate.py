#!/usr/bin/env python3
"""Aggregate Amira output directories into TSV summaries (standard library only).

Each Amira run writes `amira_results.tsv`, one row per detected AMR gene allele.
Inputs may be run directories, parent directories holding several runs, the TSV
itself, or `NAME=PATH` to choose the sample name explicitly.

Outputs
  -o/--output     per-sample summary (required)
  --long          every allele row from every sample, plus the gene family
  --matrix        sample x gene matrix (see --level, --matrix-value)
  --compare       with --amrfinder: per sample x gene family, Amira vs AMRFinderPlus

Gene naming. Amira reports tied reference names as `a/b/c` (e.g.
`blaTEM-1/blaTEM-208`). With the default `--level family` each row is collapsed
to one gene family (blaTEM-1 -> blaTEM, aac(6')-Ib-cr5 -> aac(6')-Ib) using the
FIRST name of the tie; rows whose tied names span several families are marked
`ambiguous_family` in --long. Use `--level determinant` for the verbatim names.

Copies. Amira allele counts can be inflated by read errors (one real gene split
into many low-depth alleles). `est_copies` therefore sums Amira's
"Approximate cellular copy number" per gene; compare it with `n_alleles`.
Use --min-depth / --min-identity / --min-coverage to drop weak alleles.

Examples
    amira_aggregate.py run1/ run2/ -o summary.tsv
    amira_aggregate.py A=runA/ B=runB/ -o summary.tsv --long hits.tsv --matrix m.tsv \\
        --min-depth 0.5 --amrfinder A=a.amrfinder.tsv B=b.amrfinder.tsv --compare cmp.tsv
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

RESULT_FILE = "amira_results.tsv"
NA_VALUES = {"", "N/A", "NA", "nan", "None"}
DEPTH, COPIES = "Relative mean read depth", "Approximate cellular copy number"
IDENT, COV, NAME = "Identity (%)", "Coverage (%)", "Determinant name"
BLA = r"CTX-M|TEM|OXA|SHV|KPC|CMY|NDM|VIM|IMP|EC|GES|PER|DHA|ACT|MIR|LEN|OKP"


# ---------------------------------------------------------------- helpers
def family(name):
    """Collapse an allele name to a gene family: blaTEM-1 -> blaTEM, aac(6')-Ib-cr5 -> aac(6')-Ib."""
    n = name.strip()
    n = re.sub(rf"^bla({BLA})-?.*$", r"bla\1", n)
    n = re.sub(r"-cr\d*$", "", n)
    n = re.sub(r"^(aac\(6'\)-Ib)\d+$", r"\1", n)
    return n


def families_of(determinant):
    """(chosen family, ambiguous?) for a possibly tied `a/b/c` determinant name."""
    fams = [family(x) for x in determinant.split("/") if x.strip()]
    return (fams[0] if fams else ""), len(set(fams)) > 1


def to_float(value):
    """Parse a numeric field; tied values are reported as 'a/b', take the first."""
    value = "" if value is None else str(value).strip()
    if value in NA_VALUES:
        return None
    try:
        return float(value.split("/")[0])
    except ValueError:
        return None


def mean(values):
    values = [v for v in values if v is not None]
    return sum(values) / len(values) if values else None


def fmt(value, digits=2):
    return "NA" if value is None else f"{value:.{digits}f}".rstrip("0").rstrip(".")


def split_named(arg):
    """`NAME=PATH` -> (NAME, PATH); plain PATH -> (None, PATH)."""
    if "=" in arg and not os.path.exists(arg):
        name, path = arg.split("=", 1)
        return name, path
    return None, arg


def find_result_files(paths):
    """Return [(explicit_name_or_None, sample_dir, results_tsv_or_None)]."""
    found = []
    for arg in paths:
        name, path = split_named(arg)
        if os.path.isfile(path):
            found.append((name, os.path.dirname(os.path.abspath(path)), path))
        elif not os.path.isdir(path):
            sys.stderr.write(f"WARNING: {path} does not exist, skipping.\n")
        elif os.path.exists(os.path.join(path, RESULT_FILE)):
            found.append((name, os.path.abspath(path), os.path.join(path, RESULT_FILE)))
        else:  # parent directory: search below it, do not descend into a run
            hits = 0
            for root, dirs, files in os.walk(path):
                if RESULT_FILE in files:
                    found.append((None, os.path.abspath(root), os.path.join(root, RESULT_FILE)))
                    hits += 1
                    dirs[:] = []
            if hits == 0:
                sys.stderr.write(f"WARNING: no {RESULT_FILE} found under {path}.\n")
                found.append((name, os.path.abspath(path), None))
    return found


def make_sample_names(entries):
    """Explicit NAME wins; else directory name, disambiguated with the parent, then numbered."""
    base = [n or os.path.basename(d.rstrip(os.sep)) for n, d, _ in entries]
    counts = Counter(base)
    names = []
    for (n, d, _), b in zip(entries, base):
        if counts[b] > 1 and not n:
            parent = os.path.basename(os.path.dirname(d.rstrip(os.sep)))
            b = f"{parent}_{b}" if parent else b
        names.append(b)
    seen, final = Counter(), []
    for b in names:
        seen[b] += 1
        final.append(b if seen[b] == 1 else f"{b}_{seen[b]}")
    return final


def read_results(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [r for r in reader if any((v or "").strip() for v in r.values())]
        return reader.fieldnames or [], rows


def read_amrfinder(path, scope):
    """Gene-family -> hit count from an AMRFinderPlus TSV (acquired AMR genes; scope 'core' = no plus)."""
    counts = Counter()
    with open(path, newline="") as handle:
        for r in csv.DictReader(handle, delimiter="\t"):
            if r.get("Type") != "AMR" or r.get("Subtype") != "AMR":
                continue
            if scope == "core" and r.get("Scope", "core") != "core":
                continue
            sym = r.get("Element symbol") or r.get("Gene symbol") or ""
            counts[family(sym)] += 1
    return counts


# ---------------------------------------------------------------- core
def passes(row, args):
    for col, floor in ((DEPTH, args.min_depth), (IDENT, args.min_identity), (COV, args.min_coverage)):
        if floor is None:
            continue
        v = to_float(row.get(col))
        if v is None or v < floor:
            return False
    return True


def gene_key(row, level):
    det = row.get(NAME, "")
    if level == "determinant":
        return det, False
    fam, ambiguous = families_of(det)
    return fam, ambiguous


SUMMARY_COLUMNS = [
    "sample", "status", "n_alleles", "n_alleles_filtered_out", "n_genes", "n_multi_allele_genes",
    "est_copies", "mean_identity_pct", "mean_coverage_pct", "mean_relative_depth",
    "n_flagged_alleles", "genes",
]


def summarise(sample, rows, n_dropped, level):
    per_gene = defaultdict(list)
    for r in rows:
        per_gene[gene_key(r, level)[0]].append(r)
    genes = []
    for g in sorted(per_gene):
        alleles = len(per_gene[g])
        est = sum(to_float(r.get(COPIES)) or 0 for r in per_gene[g])
        genes.append(f"{g} (alleles={alleles}, copies~{fmt(est, 1)})")
    return {
        "sample": sample,
        "status": "ok" if rows else "no_AMR_genes_detected",
        "n_alleles": len(rows),
        "n_alleles_filtered_out": n_dropped,
        "n_genes": len(per_gene),
        "n_multi_allele_genes": sum(1 for v in per_gene.values() if len(v) > 1),
        "est_copies": fmt(sum(to_float(r.get(COPIES)) or 0 for r in rows), 1),
        "mean_identity_pct": fmt(mean([to_float(r.get(IDENT)) for r in rows]), 1),
        "mean_coverage_pct": fmt(mean([to_float(r.get(COV)) for r in rows]), 1),
        "mean_relative_depth": fmt(mean([to_float(r.get(DEPTH)) for r in rows]), 2),
        "n_flagged_alleles": sum(1 for r in rows if (r.get("Comments") or "").strip()),
        "genes": "; ".join(genes),
    }


def main():
    p = argparse.ArgumentParser(description="Aggregate Amira output directories into TSV summaries.",
                                formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__.split("Examples")[1])
    p.add_argument("inputs", nargs="+", help=f"run dirs, parent dirs, {RESULT_FILE} files; optionally NAME=PATH")
    p.add_argument("-o", "--output", required=True, help="Per-sample summary TSV.")
    p.add_argument("--long", dest="long_out", help="All (kept) allele rows from all samples, with gene family.")
    p.add_argument("--matrix", dest="matrix_out", help="Sample x gene matrix TSV.")
    p.add_argument("--level", choices=["family", "determinant"], default="family",
                   help="Gene granularity for summary/matrix (default: family).")
    p.add_argument("--matrix-value", choices=["alleles", "copies"], default="alleles",
                   help="Matrix cells: allele count, or summed approximate cellular copy number (default: alleles).")
    p.add_argument("--min-depth", type=float, help="Drop alleles with relative mean read depth below this.")
    p.add_argument("--min-identity", type=float, help="Drop alleles with identity (%%) below this.")
    p.add_argument("--min-coverage", type=float, help="Drop alleles with coverage (%%) below this.")
    p.add_argument("--amrfinder", nargs="+", metavar="NAME=TSV",
                   help="AMRFinderPlus TSVs keyed by sample name (needs --compare).")
    p.add_argument("--amrfinder-scope", choices=["core", "all"], default="all",
                   help="'all' counts every AMR-type hit incl. --plus ones such as blaEC/acrF (default); 'core' keeps core-scope genes only.")
    p.add_argument("--compare", dest="compare_out", help="Per sample x gene family Amira vs AMRFinderPlus TSV.")
    args = p.parse_args()
    if bool(args.amrfinder) != bool(args.compare_out):
        p.error("--amrfinder and --compare must be used together")

    entries = find_result_files(args.inputs)
    if not entries:
        sys.exit("ERROR: no Amira results found.")
    names = make_sample_names(entries)

    summaries, long_rows, long_fields = [], [], []
    mat = defaultdict(Counter)            # sample -> gene -> alleles or copies
    amira_alleles = defaultdict(Counter)  # sample -> family -> alleles   (for --compare)
    amira_copies = defaultdict(Counter)   # sample -> family -> est copies
    genes = set()

    for (_, _, tsv), sample in zip(entries, names):
        if tsv is None:
            summaries.append({"sample": sample, "status": "missing_results_file"})
            continue
        try:
            fields, all_rows = read_results(tsv)
        except (OSError, csv.Error) as err:
            sys.stderr.write(f"WARNING: could not read {tsv}: {err}\n")
            summaries.append({"sample": sample, "status": "unreadable_results_file"})
            continue
        rows = [r for r in all_rows if passes(r, args)]
        for f in fields:
            if f not in long_fields:
                long_fields.append(f)
        summaries.append(summarise(sample, rows, len(all_rows) - len(rows), args.level))
        for r in rows:
            fam, ambiguous = families_of(r.get(NAME, ""))
            long_rows.append({"sample": sample, "gene_family": fam,
                              "ambiguous_family": "yes" if ambiguous else "no", **r})
            g = gene_key(r, args.level)[0]
            genes.add(g)
            mat[sample][g] += 1 if args.matrix_value == "alleles" else (to_float(r.get(COPIES)) or 0)
            amira_alleles[sample][fam] += 1
            amira_copies[sample][fam] += to_float(r.get(COPIES)) or 0

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=SUMMARY_COLUMNS, delimiter="\t", restval="NA")
        w.writeheader()
        w.writerows(summaries)

    if args.long_out:
        with open(args.long_out, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=["sample", "gene_family", "ambiguous_family"] + long_fields,
                               delimiter="\t", restval="NA")
            w.writeheader()
            w.writerows(long_rows)

    if args.matrix_out:
        cols = sorted(genes)
        with open(args.matrix_out, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["sample"] + cols)
            for s in names:
                w.writerow([s] + [fmt(mat[s].get(c, 0), 2) if args.matrix_value == "copies" else mat[s].get(c, 0)
                                  for c in cols])

    if args.compare_out:
        amrf = {}
        for arg in args.amrfinder:
            name, path = split_named(arg)
            if name not in names:
                sys.exit(f"ERROR: --amrfinder name {name!r} is not one of the samples: {', '.join(names)}")
            amrf[name] = read_amrfinder(path, args.amrfinder_scope)
        with open(args.compare_out, "w", newline="") as out:
            w = csv.writer(out, delimiter="\t")
            w.writerow(["sample", "gene_family", "status", "amira_alleles", "amira_est_copies", "amrfinder_hits"])
            for s in names:
                if s not in amrf:
                    continue
                for fam in sorted(set(amira_alleles[s]) | set(amrf[s])):
                    a, f = amira_alleles[s].get(fam, 0), amrf[s].get(fam, 0)
                    status = "both" if a and f else ("amira_only" if a else "amrfinder_only")
                    w.writerow([s, fam, status, a, fmt(amira_copies[s].get(fam, 0), 2), f])

    wrote = [args.output] + [x for x in (args.long_out, args.matrix_out, args.compare_out) if x]
    sys.stderr.write(f"Aggregated {len(summaries)} samples -> {', '.join(wrote)}\n")


if __name__ == "__main__":
    main()
