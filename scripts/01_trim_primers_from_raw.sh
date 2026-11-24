#!/bin/bash
# Step 1 (from raw FASTQs): remove both primers and the Nextera read-through in one pass
# using linked adapters. Run on the server, then copy only the trimmed FASTQs.
# Usage: ./01_trim_primers_from_raw.sh <raw_fastq_dir> <out_dir> <log_dir> [threads]
#
# CH16 primers: F = CGAGAAGACCCTRTGGAGCT, R = CCTNGGTCGCCCCAAC
# R1 = ^F ... revcomp(R) ; R2 = ^R ... revcomp(F)
# The 5' primer is anchored (reads start with it; no spacers in this library).
# The pair is discarded if ANY read lacks its linked primers. This also removes reads
# from the other amplicons pooled in the same index (CE16 and EU18).
set -euo pipefail
IN="$1"; OUT="$2"; LOG="$3"; THREADS="${4:-8}"; CUTADAPT=${CUTADAPT:-cutadapt}
mkdir -p "$OUT" "$LOG"
"$CUTADAPT" --version > "$LOG/cutadapt_version.txt"
for R1 in "$IN"/CH-*_R1_001.fastq.gz; do
  R2=${R1/_R1_/_R2_}; B=$(basename "$R1"); S=${B%%_R1_*}
  [ -s "$OUT/$B" ] && continue
  "$CUTADAPT" \
    -a "^CGAGAAGACCCTRTGGAGCT...GTTGGGGCGACCNAGG" \
    -A "^CCTNGGTCGCCCCAAC...AGCTCCAYAGGGTCTTCTCG" \
    -e 0.1 --discard-untrimmed --pair-filter=any -m 40 -j "$THREADS" \
    -o "$OUT/$B" -p "$OUT/$(basename "$R2")" "$R1" "$R2" > "$LOG/${S}.log" 2>&1
done
echo "done: $(ls "$OUT"/*_R1_001.fastq.gz | wc -l) samples"
