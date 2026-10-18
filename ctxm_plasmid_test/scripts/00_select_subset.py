#!/usr/bin/env python3
"""
Select the test subset from the PRJNA379782 assemblies using AllTheBacteria (ATB) precomputed tables.

Run from ctxm_plasmid_test/ after `atb fetch --tables assembly.parquet,assembly_stats.parquet,mlst.parquet,amrfinderplus.parquet`
(amrfinderplus.parquet is ~1.2 GB; do NOT let atb build its 35 GB AMR index). Requires: duckdb pandas.

Rules (see README):
  * blaCTX-M isolates: per (variant, ST) group the 2 best-assembled (fewest contigs, then highest N50) isolates;
    singleton groups are kept only for rare variants (-24, -27), plain blaCTX-M and dual-variant isolates.
  * plus 3 extra ST216 CTX-M-9 isolates (largest group), taken at ranks 1, 4, 9 of the not-yet-selected ones
    ordered by contig count, to spread assembly quality.
  * controls: best-assembled non-CTX-M isolate in each ST present in the CTX-M selection.
  * MOB-recon order: every CTX-M-15 carrier first, then controls, then the rest.

Writes: input/atb_amrfinderplus_prjna379782.tsv, input/amrfinder_atb/<SAMN>.amrfinder.tsv, input/subset.tsv,
        input/subset_biosamples.txt, input/ctxm15_isolates.tsv, input/mobrecon_order.txt
"""
import argparse
from pathlib import Path

import duckdb
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--atb-data", type=Path, default=Path.home() / ".local/share/atb/data")
ap.add_argument("--biosamples", type=Path, default=Path("input/biosamples.txt"))
ap.add_argument("--outdir", type=Path, default=Path("input"))
a = ap.parse_args()
D, O = a.atb_data, a.outdir
(O / "amrfinder_atb").mkdir(parents=True, exist_ok=True)

c = duckdb.connect()
c.sql(f"create table s as select * from read_csv('{a.biosamples}', header=false, columns={{'sid':'VARCHAR'}})")
c.sql(f"create table amr as select a.* from '{D}/amrfinderplus.parquet' a join s on a.Name = s.sid")
c.sql(f"""copy (select * from amr order by Name, "Contig id", Start) to '{O}/atb_amrfinderplus_prjna379782.tsv'
          (delimiter '\t', header)""")
c.sql(f"""create table st as select "sample" as sid, mlst_st from '{D}/mlst.parquet'""")
c.sql(f"""create table q as select sample_accession as sid, number as ncontigs, N50 from '{D}/assembly_stats.parquet'""")
c.sql("""create table ctx as select Name as sid, string_agg(distinct "Element symbol", ';' order by "Element symbol") as ctx
         from amr where "Element symbol" ilike 'blaCTX-M%' group by 1""")
m = c.sql("select ctx.*, st.mlst_st, q.ncontigs, q.N50 from ctx join q using(sid) left join st using(sid)").df()
m["mlst_st"] = m["mlst_st"].fillna("NA").astype(str)
m = m.sort_values(["ncontigs", "N50", "sid"], ascending=[True, False, True], kind="mergesort")
m["grp"] = m.ctx + "|ST" + m.mlst_st

picked = []
for g, d in m.groupby("grp", sort=True):
    k = 2 if len(d) > 1 else 1
    picked.append(d.head(k).assign(why=f"top{k} of {len(d)} in {g}"))
ct = pd.concat(picked)
rare = ct.ctx.str.contains("-24|-27|;|^blaCTX-M$")
ct = ct[~ct.why.str.contains("of 1 ") | rare]

extra = m[(m.ctx == "blaCTX-M-9") & (m.mlst_st == "216") & ~m.sid.isin(ct.sid)].sort_values(
    ["ncontigs", "sid"], kind="mergesort").iloc[[0, 3, 8]].assign(why="extra ST216 CTX-M-9 (largest group)")

nonctx = c.sql("""select q.sid, st.mlst_st, q.ncontigs, q.N50 from q left join st using(sid)
                  where q.sid in (select Name from amr) and q.sid not in (select sid from ctx)""").df()
nonctx["mlst_st"] = nonctx["mlst_st"].fillna("NA").astype(str)
ctrl = (nonctx[nonctx.mlst_st.isin(set(ct.mlst_st))].sort_values(["ncontigs", "sid"], kind="mergesort")
        .groupby("mlst_st", sort=False).head(1).assign(ctx="", role="control_no_ctxm",
                                                       why="best-assembled non-CTX-M isolate of same ST (candidate reference)"))
ct = pd.concat([ct, extra]).assign(role="ctxm")

sub = (pd.concat([ct, ctrl])[["sid", "role", "ctx", "mlst_st", "ncontigs", "N50", "why"]]
       .rename(columns={"sid": "sample", "ctx": "ctxm_variant", "mlst_st": "ST", "ncontigs": "n_contigs"})
       .sort_values(["role", "ctxm_variant", "ST", "sample"], kind="mergesort"))
sub.to_csv(O / "subset.tsv", sep="\t", index=False)
sub["sample"].to_csv(O / "subset_biosamples.txt", index=False, header=False)

c15 = m[m.ctx.str.contains("CTX-M-15")].sort_values(["ncontigs", "sid"], kind="mergesort")
c15[["sid", "ctx", "mlst_st", "ncontigs", "N50"]].rename(columns={"sid": "sample"}).to_csv(
    O / "ctxm15_isolates.tsv", sep="\t", index=False)
ctrl_ids = list(ctrl.sid)
rest = ctrl_ids + [s for s in sub["sample"] if s not in set(c15.sid) and s not in set(ctrl_ids)]
order = list(c15.sid) + rest
(O / "mobrecon_order.txt").write_text("\n".join(order) + "\n")

amr = pd.read_csv(O / "atb_amrfinderplus_prjna379782.tsv", sep="\t", low_memory=False)
for s in order:
    amr[amr.Name == s].to_csv(O / "amrfinder_atb" / f"{s}.amrfinder.tsv", sep="\t", index=False)
print(f"{len(m)} CTX-M isolates; subset {len(sub)} ({sub.role.value_counts().to_dict()}); "
      f"{len(c15)} CTX-M-15 carriers; MOB-recon order {len(order)}")
