#!/bin/bash
# HPICDA Stage 2 — download RSV Sunshine 2025 SRA runs (PRJNA1281304)
# 8 amplicon runs, ~1 GB each. prefetch -> fasterq-dump -> gzip -> verify.
set -u
export PATH="/workspace/sratoolkit.3.1.1-ubuntu64/bin:$PATH"
BASE=/workspace/HPICDA/03_functional_screens/raw/PRJNA1281304
mkdir -p "$BASE/sra" "$BASE/fastq"
RUNS="SRR34259826 SRR34259827 SRR34259828 SRR34259829 SRR34259830 SRR34259831 SRR34259832 SRR34259833"
STATUS="$BASE/download_status.tsv"
printf "run\tprefetch_bytes\tfastq_reads\tstatus\n" > "$STATUS"
for R in $RUNS; do
  echo "=== $R prefetch ==="
  prefetch --output-directory "$BASE/sra" --max-size 50G "$R" || { echo -e "$R\t0\t0\tPREFETCH_FAIL" >> "$STATUS"; continue; }
  SZ=$(stat -c%s "$BASE/sra/$R/$R.sra" 2>/dev/null || echo 0)
  echo "=== $R fasterq-dump (size=$SZ) ==="
  fasterq-dump "$BASE/sra/$R/$R.sra" --outdir "$BASE/fastq" --threads 8 --split-files || { echo -e "$R\t$SZ\t0\tDUMP_FAIL" >> "$STATUS"; continue; }
  # single-end amplicon: expect $R.fastq (or _1 if split)
  F="$BASE/fastq/$R.fastq"
  [ -f "$F" ] || F="$BASE/fastq/${R}_1.fastq"
  N=$(($(wc -l < "$F") / 4))
  gzip -1 "$F"
  echo -e "$R\t$SZ\t$N\tOK" >> "$STATUS"
  rm -rf "$BASE/sra/$R"
done
echo "ALL_DONE"
cat "$STATUS"
