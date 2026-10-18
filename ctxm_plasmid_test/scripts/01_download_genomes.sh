#!/usr/bin/env bash
# Download the short-read (Illumina) assemblies of BioProject PRJNA379782 from AllTheBacteria.
# Run from ctxm_plasmid_test/. Requires atb-cli (https://github.com/AllTheBacteria/atb-cli) and `atb fetch` done once.
set -euo pipefail
# 1. unique BioSample accessions from the SRA run table
python3 -I - <<'PY'
import pandas as pd
pd.read_csv("input/SraRunInfo_6.csv")["BioSample"].drop_duplicates().to_csv("input/biosamples.txt", index=False, header=False)
PY
# 2. resolve them to ATB assembly URLs (only the small core tables are needed; skip the 35 GB AMR index)
atb query --sample-file input/biosamples.txt --columns sample_accession,aws_url --format tsv -o input/atb_query.tsv
# 3. download
atb download --from input/atb_query.tsv --output-dir genomes --parallel 8
