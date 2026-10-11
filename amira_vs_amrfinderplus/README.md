# Amira vs AMRFinderPlus on E. coli ST131

Compares AMR gene calls from **[Amira](https://github.com/Danderson123/amira)** (directly on nanopore reads **SRR13116389**)
with **[AMRFinderPlus](https://github.com/ncbi/amr)** (on assembly **GCF_025392675.1**). Everything ran through `pixi`; the session is
logged with **[dochist](https://github.com/motroy/dochist-docs/releases/tag/v0.12.1)** v0.12.1.

> ## Read this first: the inputs are different isolates
> - **SRR13116389** — BioProject PRJNA679678 (BioSample SAMN16843782), *MinION* ST131 study of 12 strains, single-end ONT, 800k reads / 5.6 Gb.
> - **GCF_025392675.1** — BioProject PRJNA747261 (BioSample SAMN20283567), strain EcE3-19, Hospital Ramon y Cajal (Madrid), Illumina NovaSeq + SPAdes, 101 contigs.
>
> Both are *E. coli* ST131, but they are **not the same sample**. Differences in gene content below are therefore *expected biology*
> and cannot be attributed to the tools. Only genes present in both can say something about tool behaviour.

## Setup
| | |
|---|---|
| Amira | 0.11.0 (bioconda), Pandora 0.12.0-alpha.0, E. coli panRG (`Escherichia.coli.panidx.zip`), 4 cores, 1039 s |
| Reads | first **100,000 of 800,000** reads (683 Mb, ~130x, mean 6.8 kb) — subset, not the full run |
| AMRFinderPlus | 4.2.7, database 2026-08-07.1, `-O Escherichia --plus` |
| Environments | two pixi envs (`amira`, `amrfinder`): Amira's pinned flye 2.9.3 forces an old zlib that blocks AMRFinderPlus >= 4.2 |

## Results (gene-family level, `results/summary/`)
Amira: 72 allele rows → 22 gene families. AMRFinderPlus (acquired AMR genes, `Subtype=AMR`): 17 hits → 17 families.
**12 shared · 10 Amira-only · 5 AMRFinderPlus-only.**

- **Shared (12):** aac(3)-IId, aadA5, aph(3'')-Ib, aph(6)-Id, blaEC, blaOXA, blaTEM, dfrA17, mph(A), sul1, sul2, tet(A).
- **Amira-only (10):** aac(3)-IIe, aac(6')-Ib, aadA1, aph(3')-Ia, blaCMY, **blaCTX-M**, dfrA12, dfrA5, dfrA7, ere(A) — plausibly present in the reads' isolate but absent from the Madrid assembly.
- **AMRFinderPlus-only (5):** **blaKPC-31**, mrx(A) (the mph(A)-mrx(A) pair; Amira found mph(A) only), and the chromosomal efflux genes acrF, emrD, mdtM (`plus` scope; not in the NCBI acquired-gene set that Amira's panRG covers, so expected).
- AMRFinderPlus also reports gyrA/parC/parE/ptsI/uhpT point mutations and virulence/stress genes — out of Amira's scope (no point-mutation calling without `--promoter-mutations`).

### Caveat on Amira's allele/copy-number output
Amira reported many alleles per family (9 blaTEM-like, 8 blaCTX-M-like, 6 mph(A), 4 tet(A)…), each with relative depth of only ~0.4–1.2,
many at 91–99 % identity to the closest reference. For a single isolate this looks like **read-error-driven allele splitting on older R9-era nanopore data**,
not real multi-copy genes; I did not validate that (no truth assembly for the read isolate). Treat the Amira allele *counts* here as unreliable and the
*presence/absence* calls as the meaningful comparison. A fair copy-number benchmark would need a long-read assembly of the same isolate.

## Follow-up: matched-isolate comparison (same reads, assembled with Flye)
No matched pair exists publicly (the reads' isolate has no assembly or Illumina data; the assembly's isolate has no reads), so I assembled the
**same SRR13116389 reads** (first 50,000 reads, ~58x) with Flye 2.9.3 (`--nano-raw`, no polishing; 28 contigs, largest 2.0 Mb) and ran AMRFinderPlus 4.2.7 on that.
Now any difference is tool/assembly behaviour, not biology. (`results/summary_same_isolate/`, `results/amrfinderplus/SRR13116389_flye.amrfinderplus.tsv`)

**13 shared · 9 Amira-only · 5 AMRFinderPlus-only** (families, acquired AMR genes).
- Amira-only: aac(3)-IId, aph(3')-Ia, aph(6)-Id, blaCMY, dfrA5/7/12, ere(A), sul2.
- AMRFinderPlus-only: catB3, catB8, fosA, dfrA (unresolved subtype), mdtM (`plus`). Amira's panRG may not contain these, or it filtered them; I did not check which.
- **18 of the 32 AMR-type hits on the Flye assembly are `PARTIALX`** (50–85 % coverage of the reference), consistent with uncorrected nanopore indels breaking
  gene frames in an unpolished assembly; this probably explains some of the Amira-only genes. I did not test this (e.g. by polishing with medaka/racon), so it is a hypothesis.
- Flye split many genes across multiple contigs/copies (e.g. sul1 x3, aac(6')-Ib x3), so AMRFinderPlus copy counts here are also not trustworthy as truth.
- Caveat still applies: no ground-truth assembly, so neither tool is "right"; Amira's high per-family allele counts remain unvalidated.

## Layout
- `pixi.toml`, `pixi.lock` — environments
- `scripts/run_all.sh` — clean end-to-end reproduction; `scripts/compare.py` — comparison
- `results/amrfinderplus/` — AMRFinderPlus table + detected gene sequences
- `results/amira_curated/` — `amira_results.tsv`, final polished allele FASTA, coverage plot, run log (the 2.7 GB raw Amira output dir is git-ignored)
- `results/summary/` — cross-isolate comparison table; `results/summary_same_isolate/` — Amira vs Flye+AMRFinderPlus on the same reads; `results/flye_curated/` — Flye assembly
- `amira-vs-amrfinderplus.dochist.json` — saved dochist session (16 commands, includes the failed/superseded steps); `FAIR-report.md` — dochist FAIR report. Reload with `dochist load amira-vs-amrfinderplus.dochist.json`.

## Notes on the run
- `figshare.com/ndownloader/...` returns a bot challenge (HTTP 202, empty body); the same file downloads normally from `ndownloader.figshare.com`.
- The `pypi` package `amira-amr==0.11.0` did not resolve; used the bioconda `amira` package.
- The first AMRFinderPlus run used 3.12.8 / DB 2024-07-22.1 (solver default); it was redone with 4.2.7 / DB 2026-08-07.1. Both attempts are in the dochist log.
