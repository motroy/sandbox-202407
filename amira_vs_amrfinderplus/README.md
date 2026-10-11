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

## Matched-isolate benchmark: SAMN08768108 (*E. coli* O103:H2 FWSEC0007), myloasm + Polypolish + pypolca
A true matched set, with a **complete genome for truth** (GCF_005037795.1: chromosome + 1 plasmid): ONT MinION SRR7896269 (100,004 reads, 312 Mb, mean 3.1 kb, ~57x)
and MiSeq 2x SRR7896254 (1.28 M pairs, 567 Mb, ~100x). Pipeline: `scripts/fwsec0007.sh`; summary: `scripts/fwsec_summary.py`; outputs in `results/FWSEC0007/`.
Tools: myloasm 0.7.0 (defaults) → Polypolish 0.7.1 (bwa `-a`, `filter`, `polish`) → pypolca 0.5.0 (`--careful`); AMRFinderPlus 4.2.7; Amira 0.11.0 (323 s).

**The truth genome has no acquired AMR genes.** AMRFinderPlus on it reports only chromosomal acrF, mdtM and blaEC (`plus` scope), plus stress and a stx1a virulence operon.

| assembly stage | contigs | total | error vs truth (aligned blocks) |
|---|---|---|---|
| myloasm | 272 | 7.94 Mb | 1.77 % (≈QV 17.5) |
| + Polypolish | 272 | 7.96 Mb | 0.65 % (≈QV 21.9) |
| + pypolca | 272 | 7.97 Mb | 0.53 % (≈QV 22.7) |

| AMRFinderPlus / Amira call | truth | myloasm | +polypolish | +pypolca | Amira |
|---|---|---|---|---|---|
| acrF, mdtM (chromosomal efflux) | ✓ | ✗ | ✓ | ✓ | – (not an acquired gene) |
| blaEC | ✓ | ✗ | ✓ | ✓ | ✓ (blaEC-18, 1 allele) |
| spurious frameshift/“POINT_DISRUPT” calls (acrR, marR, nfsA, ompC) | 0 | 3 | 2 | 2 | 0 |

- **Raw myloasm assembly loses every gene**: indels break the reading frames (AMRFinderPlus reports FRAMESHIFT and no BLASTX hits for acrF/mdtM/blaEC). **Polypolish alone restores them**; pypolca then removes one spurious call (marR) but leaves **2 false frameshift “resistance” calls** (nfsA, ompC) that the truth genome does not have.
- **Amira recovered the one relevant acquired-style gene (blaEC) with a single allele at relative depth 1.0 and no false positives** from the unassembled reads. This is much cleaner than the ST131 run above, which suggests that run's many-allele output was driven by that sample's multi-copy-rich/noisy content rather than a general Amira problem — still not proven.
- **The myloasm assembly is not good here**: 272 contigs / 7.9 Mb for a 5.5 Mb genome (N50 66 kb, none ≥ 1 Mb; no circular chromosome). myloasm targets R10 reads (>97 % accuracy); these are R9-era reads (≈90 % mean identity estimated from base qualities), and the reads are short (N50 5.4 kb). Residual error after polishing is still 0.53 % (≈QV 23) — far from a finished genome.
- Trials: `--quality-value-cutoff 85` (keeps 93 % of reads instead of 58 %) was **worse** (359 contigs, 9.5 Mb, N50 43 kb), so defaults were kept.
- Method caveat: error rate is `NM / aligned block length` from `minimap2 -cx asm5` against truth. Because the assembly is redundant, aligned bases exceed what a single copy would give, so treat it as a rough comparison between stages, not a formal QV.
- This isolate carries essentially no acquired resistance, so it tests false positives and recovery of intrinsic genes, not detection of real resistance genes; a sample with plasmid-borne AMR genes, ideally R10 ONT, would be a more informative test of myloasm and Amira.

### Same sample, assembled with **dragonflye** instead of myloasm
Dragonflye 1.2.1 (Flye 2.9.x `--nano-raw` + 1 Racon round + dnaapler reorientation; `--gsize 5.5m`, other defaults, no Medaka), then the same Polypolish → pypolca steps
(`scripts/fwsec0007.sh dragonflye | polish_df | amrfinder_df`; 9.5 min assembly). Accuracy is against the complete genome, same method and caveats as above
(`results/FWSEC0007/summary/`).

| stage | contigs | total | truth covered | error vs truth | ≈QV |
|---|---|---|---|---|---|
| myloasm | 272 | 7.94 Mb | 86.0 % | 1.77 % | 17.5 |
| myloasm + Polypolish + pypolca | 272 | 7.97 Mb | 90.2 % | 0.53 % | 22.7 |
| dragonflye (Flye) | 12 | 5.43 Mb | 93.8 % | 1.60 % | 17.9 |
| dragonflye (Flye + Racon) | 4 | 5.45 Mb | 100 % | 0.79 % | 21.0 |
| dragonflye + Polypolish | 4 | 5.48 Mb | 100 % | 0.037 % | 34.3 |
| **dragonflye + Polypolish + pypolca** | **4** | **5.48 Mb** | **100 %** | **0.030 %** | **35.2** |

- Dragonflye's final assembly is 4 contigs (chromosome split in two, 3.76 + 1.62 Mb; one circular 72.8 kb plasmid-sized contig; one 1.5 kb low-coverage contig), versus 272 for myloasm on the same R9-era reads.
- **AMR calls from the final dragonflye assembly match the truth genome exactly** (acrF, mdtM, blaEC) with **no spurious frameshift calls**. The raw Flye and Racon stages show false “POINT_DISRUPT” calls
  (acrR, marR, cirA) and miss all three genes; **Polypolish is what fixes it** (Racon 0.79 % → Polypolish 0.037 % error); pypolca adds a small further gain (0.037 → 0.030 %). Amira's single blaEC call is consistent with this.
- One isolate, one run per tool, ~57x of older ONT chemistry: this is a worked example, not a benchmark of the assemblers in general. myloasm is built for R10 reads, so this comparison is unfavourable to it.

## `scripts/amira_aggregate.py`: summarise many Amira runs
Stdlib-only aggregator (extended from a user-supplied script). Collapses tied names (`blaTEM-1/blaTEM-208`) to gene families, reports `est_copies`
(sum of Amira's approximate cellular copy number) beside the raw allele count, filters weak alleles (`--min-depth/--min-identity/--min-coverage`),
and compares against AMRFinderPlus (`--amrfinder NAME=TSV --compare out.tsv`, `both`/`amira_only`/`amrfinder_only`).
```
python -I scripts/amira_aggregate.py ST131=results/amira/SRR13116389 FWSEC0007=results/FWSEC0007/amira \
  -o summary.tsv --long hits.tsv --matrix matrix.tsv --matrix-value copies \
  --amrfinder ST131=results/amrfinderplus/SRR13116389_flye.amrfinderplus.tsv FWSEC0007=results/FWSEC0007/dragonflye/amrfinder/polished.tsv --compare compare.tsv
```
