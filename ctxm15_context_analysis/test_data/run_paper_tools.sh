#!/usr/bin/env bash
# Run the paper's tools on the test subset (run from ctxm15_context_analysis/):
#   ABRicate (ResFinder, ISfinder, FimTyper, EcOH; 80% identity/coverage) and Platon, in test_data/run_paper/.
# Needs abricate (1.0.1 in the paper) with custom `isfinder` and `fimtyper` databases (see README), platon + PLATON_DB.
set -euo pipefail
R=test_data/run_paper; mkdir -p $R/abricate $R/platon
THREADS=${THREADS:-1}; JOBS=${JOBS:-4}
samples=$(cut -f2 test_data/subset.tsv | tail -n +2)
for s in $samples; do
  for db in resfinder isfinder fimtyper ecoh; do
    [ -s $R/abricate/$s.$db.tsv ] || abricate --db $db --minid 80 --mincov 80 --quiet test_data/genomes/$s.fa > $R/abricate/$s.$db.tsv
  done
done
if [ -n "${PLATON_DB:-}" ]; then
  echo "$samples" | xargs -P "$JOBS" -I{} bash -c '[ -s '$R'/platon/{}/{}.tsv ] || platon --db "$PLATON_DB" --threads '"$THREADS"' --output '$R'/platon/{} test_data/genomes/{}.fa >/dev/null 2>&1'
else
  echo "[skip] Platon: set PLATON_DB" >&2
fi
