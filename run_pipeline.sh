#!/bin/bash
# Run the full pipeline. Usage: ./run_pipeline.sh <raw_fastq_dir> <reference_fasta> [threads]
set -euo pipefail
RAW="$1"; REF="$2"; THREADS="${3:-8}"
mkdir -p work logs
./scripts/01_trim_primers_from_raw.sh "$RAW" work/trimmed logs/cutadapt "$THREADS"
Rscript scripts/02_dada2.R work/trimmed work/filtered . "$REF" "$THREADS"
Rscript scripts/02b_gulf_checklist_filter.R results/CH_taxonomy_minBoot80.csv metadata/gulf_of_california_fish_checklist.csv results/CH_taxonomy_minBoot80_curated.csv
Rscript scripts/03_diet_analysis.R . metadata/sample_metadata.csv metadata/gulf_of_california_fish_checklist.csv NA 999
Rscript scripts/04_compare_read_thresholds.R analysis metadata/sample_metadata.csv
