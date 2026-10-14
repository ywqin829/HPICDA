#!/bin/bash
# HPICDA Stage 2 — RSV PRJNA1281304: mageck count all 8 runs (combined GeCKOv2 A+B)
set -u
export PATH="/workspace/.conda/envs/crispr/bin:$PATH"
BASE=/workspace/HPICDA/03_functional_screens/raw/PRJNA1281304
cd "$BASE"
# sample name mapping
declare -A NAME=( [SRR34259833]=CTRL_A [SRR34259832]=CTRL_B \
  [SRR34259831]=A_GFP_LOW [SRR34259830]=A_GFP_NEG \
  [SRR34259829]=B_GFP_LOW [SRR34259828]=B_GFP_NEG \
  [SRR34259827]=viability_A [SRR34259826]=viability_B )
for R in SRR34259826 SRR34259827 SRR34259828 SRR34259829 SRR34259830 SRR34259831 SRR34259832 SRR34259833; do
  N=${NAME[$R]}
  echo "=== counting $R ($N) ==="
  mageck count -l Human_GeCKOv2_Library_combine.csv -n "count_$N" \
    --fastq "fastq/$R.fastq.gz" --trim-5 0 --sgrna-len 20 2>&1 | grep -iE "Mapped:|Total:|error"
done
echo "COUNT_ALL_DONE"
