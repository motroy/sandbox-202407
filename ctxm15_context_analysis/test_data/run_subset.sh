#!/usr/bin/env bash
# Run the pipeline on the Zenodo test subset. Run from ctxm15_context_analysis/.
# Needs amrfinder, mob_recon, mlst, ectyper on PATH, genomes in test_data/genomes (download_genomes.py)
# and the K-12 MG1655 GBFF (NC_000913.3) at test_data/k12_MG1655.gbff.
set -euo pipefail
R=test_data/run; mkdir -p $R/{amrfinder,mobrecon,typing}
THREADS=${THREADS:-2}; JOBS=${JOBS:-4}
cut -f2 test_data/subset.tsv | tail -n +2 | xargs -P "$JOBS" -I{} bash -c '
  s={}; fa=test_data/genomes/$s.fa; cd '"$R"'
  [ -s amrfinder/$s.amrfinder.tsv ] || amrfinder -n ../../$fa --plus --threads '"$THREADS"' -o amrfinder/$s.amrfinder.tsv
  [ -s mobrecon/$s/contig_report.txt ] || mob_recon -i ../../$fa -o mobrecon/$s -n '"$THREADS"' --force
  [ -s typing/$s.mlst.tsv ] || mlst ../../$fa > typing/$s.mlst.tsv
  [ -s typing/$s.ectyper.tsv ] || { ectyper -i ../../$fa -o typing/$s.ectyper --cores '"$THREADS"' && cp typing/$s.ectyper/output.tsv typing/$s.ectyper.tsv; }
'
python3 - <<'PY'
import csv
rows = list(csv.DictReader(open("test_data/subset.tsv"), delimiter="\t"))
with open("test_data/run/samples.tsv", "w") as f:
    f.write("sample\tfasta\tamrfinder\tmobrecon_dir\tst\tgbff\tmlst\tclermont\tectyper\tfimtyper\n")
    for r in rows:
        s = r["Strain"]; st = r["Sequence Type"]
        f.write(f"{s}\ttest_data/genomes/{s}.fa\ttest_data/run/amrfinder/{s}.amrfinder.tsv\ttest_data/run/mobrecon/{s}\t"
                f"{'ST' + st if st != 'ND' else 'ND'}\ttest_data/genomes/{s}.gbff\ttest_data/run/typing/{s}.mlst.tsv\t\t"
                f"test_data/run/typing/{s}.ectyper.tsv\t\n")
PY
python3 scripts/02_context_analysis.py --config test_data/config.test.toml
[ -s test_data/k12_MG1655.gbff ] && python3 scripts/03_insertion_sites.py test_data/config.test.toml || echo "[skip] 03: no K-12 GBFF"
python3 scripts/05_typing.py test_data/config.test.toml
python3 test_data/compare_to_paper.py test_data/subset.tsv test_data/run/results
