#!/bin/bash
# Step 1: remove the 3' primer (reverse complement) plus the Nextera adapter read-through.
# Input: chordata/for_dada2 (5' primers already removed with cutadapt -g/-G --discard-untrimmed).
# R1 3' end: reverse complement of CH_R (CCTNGGTCGCCCCAAC) = GTTGGGGCGACCNAGG
# R2 3' end: reverse complement of CH_F (CGAGAAGACCCTRTGGAGCT) = AGCTCCAYAGGGTCTTCTCG
# Pairs are discarded if either read lacks its 3' primer (the ~77 bp insert is always read through).
set -euo pipefail
IN="$1"; OUT="$2"; LOG="$3"; CUTADAPT=${CUTADAPT:-cutadapt}
mkdir -p "$OUT" "$LOG"
"$CUTADAPT" --version > "$LOG/cutadapt_version.txt"
for R1 in "$IN"/CH-*_R1_001.fastq.gz; do
  R2=${R1/_R1_/_R2_}; B=$(basename "$R1"); S=${B%%_R1_*}
  [ -s "$OUT/$B" ] && continue
  "$CUTADAPT" -a GTTGGGGCGACCNAGG -A AGCTCCAYAGGGTCTTCTCG \
    -O 10 -e 0.1 --discard-untrimmed --pair-filter=any -m 40 -j 8 \
    -o "$OUT/$B" -p "$OUT/$(basename "$R2")" "$R1" "$R2" > "$LOG/${S}.log" 2>&1
done
echo done
