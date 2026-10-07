"""Shared helpers: MetaPhlAn / HUMAnN parsing, condition mapping, lineage handling, logging."""
import glob
import gzip
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", message="This pattern is interpreted as a regular expression")

RANK_PREFIX = {"k": "Kingdom", "p": "Phylum", "c": "Class", "o": "Order", "f": "Family", "g": "Genus", "s": "Species"}

# Suffixes stripped from file names / column headers to obtain sample ids
DEFAULT_SAMPLE_REGEX = (
    r"((_profile|_metaphlan\d*|_bugs_list|_pathabundance(_cpm|_relab)?|_Abundance(-RPKs|-CPM|-RELAB)?)"
    r"|(\.(txt|tsv|csv)(\.gz)?))+$"
)


def setup_log(snakemake):
    """Send stdout/stderr of a script: job to its log file."""
    if snakemake.log:
        sys.stderr = sys.stdout = open(snakemake.log[0], "w", buffering=1)


def open_any(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def strip_sample_id(name, regex=None):
    return re.sub(regex or DEFAULT_SAMPLE_REGEX, "", os.path.basename(str(name)))


def expand_paths(pattern):
    """A glob, a directory-less file path, or a list of either -> sorted list of existing files."""
    patterns = [pattern] if isinstance(pattern, str) else list(pattern)
    files = []
    for p in patterns:
        hits = sorted(glob.glob(p))
        files.extend(hits)
    return files


# ---------------------------------------------------------------------------
# Lineage helpers (features are MetaPhlAn clade names: k__A|p__B|...|s__Species)
# ---------------------------------------------------------------------------
def species_key(clade):
    return clade.split("|")[-1]


def parse_lineage(clade):
    """'k__Bacteria|p__Firmicutes|...' -> {'Kingdom': 'Bacteria', 'Phylum': 'Firmicutes', ...}"""
    out = {}
    for part in clade.split("|"):
        if len(part) > 3 and part[1:3] == "__" and part[0] in RANK_PREFIX:
            out[RANK_PREFIX[part[0]]] = part[3:]
    return out


# ---------------------------------------------------------------------------
# MetaPhlAn (v3/v4) readers -> DataFrame features x samples, species level, percent
# ---------------------------------------------------------------------------
def _read_metaphlan_file(path, sample_regex=None):
    header, rows = None, []
    with open_any(path) as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if not line:
                continue
            if line.startswith("#"):
                if line.startswith("#clade_name"):
                    header = line[1:].split("\t")
                continue
            if line.startswith("clade_name"):
                header = line.split("\t")
                continue
            rows.append(line.split("\t"))
    if not rows:
        raise ValueError(f"No data rows in MetaPhlAn file {path}")

    per_sample = header is None or "relative_abundance" in header or "NCBI_tax_id" in header
    if per_sample:
        if header and "relative_abundance" in header:
            idx = header.index("relative_abundance")
        else:
            idx = 2 if len(rows[0]) >= 3 else 1
        sid = strip_sample_id(path, sample_regex)
        data = {r[0]: float(r[idx]) for r in rows if len(r) > idx and r[idx] != ""}
        return pd.DataFrame({sid: pd.Series(data)})
    # merged table (merge_metaphlan_tables.py): clade_name + one column per sample
    cols = [strip_sample_id(c, sample_regex) for c in header[1:]]
    df = pd.DataFrame([r[1:1 + len(cols)] for r in rows], index=[r[0] for r in rows], columns=cols)
    return df.apply(pd.to_numeric, errors="coerce")


def read_metaphlan(pattern, sample_regex=None):
    """Read per-sample profiles and/or merged tables; keep species-level clades (s__, no t__ SGB leaves)."""
    files = expand_paths(pattern)
    if not files:
        raise FileNotFoundError(f"No MetaPhlAn files match {pattern!r}")
    frames = []
    for f in files:
        df = _read_metaphlan_file(f, sample_regex)
        keep = [i for i in df.index if species_key(i).startswith("s__")]
        frames.append(df.loc[keep])
    out = pd.concat(frames, axis=1, sort=False).fillna(0.0)
    out = out.groupby(level=0).sum()
    if out.columns.duplicated().any():
        dup = out.columns[out.columns.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Duplicate sample ids in MetaPhlAn input: {dup}")
    return out


# ---------------------------------------------------------------------------
# HUMAnN pathway abundance readers -> DataFrame pathways x samples (unstratified)
# ---------------------------------------------------------------------------
def read_humann(pattern, sample_regex=None):
    files = expand_paths(pattern)
    if not files:
        raise FileNotFoundError(f"No HUMAnN files match {pattern!r}")
    frames = []
    for f in files:
        df = pd.read_csv(f, sep="\t", header=0, index_col=0, comment=None, low_memory=False)
        df.columns = [strip_sample_id(c, sample_regex) for c in df.columns]
        df.index = df.index.astype(str)
        df = df.loc[~df.index.str.contains(r"\|", regex=True)]  # unstratified rows only
        frames.append(df.apply(pd.to_numeric, errors="coerce"))
    out = pd.concat(frames, axis=1, sort=False).fillna(0.0)
    out = out.groupby(level=0).sum()
    if out.columns.duplicated().any():
        dup = out.columns[out.columns.duplicated()].unique().tolist()[:5]
        raise ValueError(f"Duplicate sample ids in HUMAnN input: {dup}")
    return out


def read_metadata(path, sample_id_column):
    sep = "," if str(path).endswith((".csv", ".csv.gz")) else "\t"
    meta = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=True)
    if sample_id_column not in meta.columns:
        raise KeyError(f"sample_id_column {sample_id_column!r} not in metadata columns {list(meta.columns)[:15]}")
    meta = meta.rename(columns={sample_id_column: "sample_id"})
    if meta["sample_id"].duplicated().any():
        raise ValueError(f"Duplicated sample ids in {path}")
    return meta


# ---------------------------------------------------------------------------
# Condition mapping
# ---------------------------------------------------------------------------
def map_condition(values, spec, labels):
    """Map a metadata column to control/case labels using *_values lists and/or *_regex patterns.
    Values matching neither (or both) become NaN."""
    s = values.astype("string")

    def match(vals_key, rx_key):
        m = pd.Series(False, index=s.index)
        if spec.get(vals_key):
            wanted = {str(v).lower() for v in spec[vals_key]}
            m |= s.str.lower().isin(wanted).fillna(False)
        if spec.get(rx_key):
            m |= s.str.contains(spec[rx_key], flags=re.IGNORECASE, regex=True).fillna(False)
        return m

    is_ctl, is_case = match("control_values", "control_regex"), match("case_values", "case_regex")
    out = pd.Series(np.nan, index=s.index, dtype=object)
    out[is_ctl & ~is_case] = labels["control"]
    out[is_case & ~is_ctl] = labels["case"]
    return out


def read_matrix(path):
    return pd.read_csv(path, sep="\t", index_col=0)


def write_matrix(df, path):
    df.index.name = "feature"
    df.to_csv(path, sep="\t")
