#!/usr/bin/env python3
"""Download FASTA + GBFF for the subset from the NCBI Datasets API (no extra tools needed).

Usage: download_genomes.py subset.tsv outdir   -> outdir/<acc>.fa and outdir/<acc>.gbff
"""
import csv, io, sys, time, urllib.request, zipfile
from pathlib import Path

URL = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
       "?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GBFF&filename={acc}.zip")
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    acc = r["Strain"]
    if (out / f"{acc}.fa").exists() and (out / f"{acc}.gbff").exists():
        continue
    for attempt in range(4):
        try:
            z = zipfile.ZipFile(io.BytesIO(urllib.request.urlopen(URL.format(acc=acc), timeout=300).read()))
            break
        except Exception as e:
            print(f"[retry {attempt + 1}] {acc}: {e}", file=sys.stderr); time.sleep(2 ** (attempt + 1))
    else:
        print(f"[FAIL] {acc}", file=sys.stderr); continue
    for n in z.namelist():
        if n.endswith(".fna"): (out / f"{acc}.fa").write_bytes(z.read(n))
        elif n.endswith(".gbff"): (out / f"{acc}.gbff").write_bytes(z.read(n))
    print(f"{acc} ok", file=sys.stderr)
