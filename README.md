# sealion-diet-metabarcoding

Bioinformatic and statistical pipeline for the fecal eDNA metabarcoding diet study of the California sea lion (*Zalophus californianus*) at Los Islotes, Gulf of California:

> Mejía-Guarnizo B. N. et al. *Fecal eDNA metabarcoding recovers documented temporal dietary patterns in California sea lions (Zalophus californianus).* (in review)

The pipeline processes 16S rRNA mitochondrial amplicons (Chord_16S primers; Deagle et al. 2009), builds an amplicon sequence variant (ASV) table with DADA2, assigns taxonomy against a GenBank-derived chordate 16S reference, and computes prey frequency of occurrence (FO), PERMANOVA, PERMDISP, balanced subsampling, a rarefaction check and a richness–depth correlation.

Raw FASTQ files and the reference database are not stored in this repository (see [Data](#data)).

## Amplicon

| | Sequence (5'→3') |
|---|---|
| Forward (CH16_F) | `CGAGAAGACCCTRTGGAGCT` |
| Reverse (CH16_R) | `CCTNGGTCGCCCCAAC` |
| Insert without primers | ~77 bp (main peak 76–79 bp; secondary peak 94 bp, *Ophichthus*) |

Libraries were sequenced 2 × 150 bp on an Illumina HiSeq X (patterned flow cell, binned quality scores) with dual Nextera XT indices. Each library pools three amplicons (chordate CH16, cephalopod CE16 and eukaryote EU18); reads are assigned to the chordate amplicon by its primers. Because the insert is much shorter than the reads, every read runs through the insert into the opposite primer and the Nextera adapter, so **both primers and the adapter read-through must be removed at both ends**.

## Pipeline

| Step | Script | What it does |
|---|---|---|
| 1 | `scripts/01_trim_primers_from_raw.sh` | cutadapt with linked adapters: anchored 5' primer + reverse complement of the opposite primer at the 3' end, on both reads (`-e 0.1 --discard-untrimmed --pair-filter=any -m 40`). This also separates the chordate amplicon from the CE16/EU18 amplicons in the same library. |
| 2 | `scripts/02_dada2.R` | `filterAndTrim` (no truncation, `maxEE = c(2,2)`, `truncQ = 2`, `maxN = 0`, `minLen = 50`, `rm.phix = TRUE`); `learnErrors` with a loess fit forced to be monotonic (for binned quality scores); per-sample `dada`; `mergePairs` (`minOverlap = 20`, `maxMismatch = 0`, `trimOverhang = TRUE`); `removeBimeraDenovo` (consensus); insert-length filter 60–100 bp; SHA-1 hashes; `assignTaxonomy` (`tryRC = TRUE`, bootstraps stored), exported at `minBoot = 80` (main) and `minBoot = 50`. |
| 3 | `scripts/03_diet_analysis.R` | Taxonomic curation, FO, PERMANOVA/PERMDISP, balanced subsampling, rarefaction to 4,000 prey reads, Spearman correlation between richness and depth. |
| 4 | `scripts/04_compare_read_thresholds.R` | Sensitivity check: presence with ≥ 1 read vs ≥ 10 reads (minBoot 80). |
| extra | `scripts/extra/01b_trim_3prime_from_5prime_trimmed.sh` | Alternative step 1 for reads whose 5' primers were already removed (used to cross-check step 1). |

### Curation and detection rules (step 3)

- Species-level assignments within Actinopteri that are not in the Gulf of California fish checklist (`metadata/gulf_of_california_fish_checklist.csv`; Robertson et al. 2024) are reported at genus level (`Genus sp`). ASVs without a genus assignment are excluded from prey tables.
- Non-prey taxa are removed: the predator (*Zalophus*), humans, domestic vertebrates (*Bos*, *Gallus*, *Sus*, *Canis*) and freshwater taxa (*Oreochromis*).
- Two detection thresholds per sample:
  - **original**: ≥ 1 read;
  - **strict** (main result): ≥ 0.1% of the sample's prey reads **and** ≥ 10 reads.
- Multivariate analyses use presence/absence, Jaccard distance and 999 permutations, with and without removing prey present in < 5% of the 96 samples. There are two one-factor comparisons: interannual (post-reproductive 2021, n = 8, vs post-reproductive 2022, n = 23) and seasonal within 2022 (reproductive, n = 65, vs post-reproductive, n = 23). PERMDISP uses `betadisper` + `permutest`. Balanced subsampling runs 999 iterations, with the larger group subsampled to the size of the smaller one.

## Usage

```bash
# 1. primer/adapter trimming (raw FASTQs named CH-<sample>_S<n>_L001_R[12]_001.fastq.gz)
CUTADAPT=cutadapt ./scripts/01_trim_primers_from_raw.sh raw_fastqs/ work/trimmed/ logs/cutadapt/ 8

# 2. DADA2 + taxonomy (writes results/ and logs/ inside <out_dir>)
Rscript scripts/02_dada2.R work/trimmed/ work/filtered/ . reference/16S_chor_taxonomy.fasta 8

# 3. diet analysis (4th argument: folder with the 2024 pipeline outputs for comparison, or NA)
Rscript scripts/03_diet_analysis.R . metadata/sample_metadata.csv metadata/gulf_of_california_fish_checklist.csv NA 999

# 4. >=1 vs >=10 reads
Rscript scripts/04_compare_read_thresholds.R analysis/ metadata/sample_metadata.csv
```

`run_pipeline.sh` chains the four steps.

## Repository contents

| Folder | Contents |
|---|---|
| `scripts/` | Pipeline scripts (steps 1–4). |
| `metadata/` | Sample metadata (`sample_metadata.csv`: year, reproductive season, site), Gulf of California fish checklist, primer table. |
| `results/` | ASV table (`CH_ASV_table.csv`, long format: sample, hash, reads), hash → sequence key (`.csv`, `.fasta`), taxonomy at minBoot 80 and 50, bootstrap values. |
| `analysis/` | Prey read matrices, FO tables, PERMANOVA/PERMDISP, subsampling, Spearman, Figure 8 (FO > 5%) and NMDS. Files named `original_*` re-analyse the 2024 pipeline output for comparison. |
| `logs/` | Per-sample read tracking, parameters, quality profiles, error models, software versions. |

## Main results (minBoot 80, strict threshold)

- 163.3 M chordate read pairs → 143.1 M reads in 2,693 ASVs (87.6% retained); 106.6 M prey reads; 74 prey taxa.
- Interannual: PERMANOVA F = 3.72, R² = 0.114, p = 0.001; PERMDISP p = 0.73.
- Seasonal: PERMANOVA F = 3.98, R² = 0.044, p = 0.001; PERMDISP p = 0.26.

## Software

| Tool | Version |
|---|---|
| cutadapt | 4.4 |
| R | 4.3.3 |
| dada2 | 1.28.0 |
| ShortRead | 1.58.0 |
| vegan | 2.6-10 |
| dplyr / tidyr / readr / tibble | 1.1.4 / 1.3.1 / 2.1.5 / 3.2.1 |
| digest | 0.6.37 |

Full session information: `logs/sessionInfo.txt` and `analysis/sessionInfo.txt`.

## Data

- **Raw reads:** [SRA/ENA accession to be added].
- **Reference database:** GenBank chordate mitochondrial sequences, in-silico PCR with the Chord_16S primers (`insect::virtualPCR`), NCBI taxonomy (`taxize`). 142,689 sequences, 21,686 species. [Zenodo DOI to be added.]

## License

Code: MIT (see `LICENSE`).
