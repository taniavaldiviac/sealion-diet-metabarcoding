<h1 align="center">sealion-diet-metabarcoding</h1>

<p align="center">
  <b>Fecal eDNA metabarcoding pipeline for the diet of the California sea lion (<i>Zalophus californianus</i>) in the Gulf of California</b>
</p>

<p align="center">
  <a href="LICENSE"><img alt="License: MIT" src="https://img.shields.io/badge/license-MIT-green.svg"></a>
  <img alt="R 4.3.3" src="https://img.shields.io/badge/R-4.3.3-276DC3?logo=r&logoColor=white">
  <img alt="DADA2 1.28.0" src="https://img.shields.io/badge/DADA2-1.28.0-orange">
  <img alt="cutadapt 4.4" src="https://img.shields.io/badge/cutadapt-4.4-blue">
  <img alt="vegan 2.6-10" src="https://img.shields.io/badge/vegan-2.6--10-brightgreen">
  <img alt="Marker: 16S rRNA" src="https://img.shields.io/badge/marker-16S%20rRNA%20(Chord__16S)-lightgrey">
  <img alt="Status: in review" src="https://img.shields.io/badge/manuscript-in%20review-yellow">
</p>

---

## Contents

- [Overview](#overview)
- [Key results](#key-results)
- [Pipeline at a glance](#pipeline-at-a-glance)
- [Amplicon and sequencing](#amplicon-and-sequencing)
- [Reference database](#reference-database)
- [Gulf of California fish filter](#gulf-of-california-fish-filter)
- [Quick start](#quick-start)
- [Step-by-step parameters](#step-by-step-parameters)
- [Curation and detection rules](#curation-and-detection-rules)
- [Repository structure](#repository-structure)
- [Software](#software)
- [Data availability](#data-availability)
- [How to cite](#how-to-cite)
- [References](#references)
- [License](#license)

## Overview

This repository contains the complete bioinformatic and statistical workflow behind:

> Mejía-Guarnizo B. N. *et al.* **Fecal eDNA metabarcoding recovers documented temporal dietary patterns in California sea lions (*Zalophus californianus*).** *In review.*

We analysed **96 scats** collected at the Los Islotes colony (Gulf of California) during three periods: post-reproductive 2021 (n = 8), post-reproductive 2022 (n = 23) and reproductive 2022 (n = 65). Prey DNA was amplified with the vertebrate **Chord_16S** primers (mitochondrial 16S rRNA, ~77 bp insert). The pipeline:

1. removes both primers and the adapter read-through at **both ends** of every read;
2. infers ASVs with **DADA2** and assigns taxonomy against a GenBank chordate 16S reference (minBoot 80);
3. applies a **sample-level detection threshold** that guards against cross-contamination and index hopping;
4. computes **frequency of occurrence**, **PERMANOVA**, **PERMDISP**, balanced subsampling, a rarefaction check and a richness–depth correlation.

## Key results

Main analysis: minBoot 80, a prey counted as present when it makes up ≥ 0.1% of the sample's prey reads and has ≥ 10 reads.

| | Value |
|---|---|
| Chordate read pairs → reads in final ASV table | 163.3 M → 143.1 M (87.6%) |
| ASVs | 2,693 |
| Prey reads / prey taxa | 106.6 M / 74 taxa (33 species, 41 genera) |
| Interannual (2021 vs 2022, post-reproductive) | PERMANOVA F = 3.72, R² = 0.114, p = 0.001 · PERMDISP p = 0.73 |
| Seasonal (2022, reproductive vs post-reproductive) | PERMANOVA F = 3.98, R² = 0.044, p = 0.001 · PERMDISP p = 0.26 |

<p align="center">
  <img src="analysis/NMDS_main.png" alt="NMDS of prey composition (Jaccard)" width="48%">
  <img src="analysis/Fig8_FO_main.png" alt="Frequency of occurrence of prey with FO > 5%" width="48%">
</p>
<p align="center"><sub>Left: NMDS of prey presence/absence (Jaccard) by period. Right: frequency of occurrence of prey detected in more than 5% of scats.</sub></p>

## Pipeline at a glance

```mermaid
flowchart LR
    R0[00 · Reference database<br/>GenBank → in-silico PCR<br/>→ NCBI taxonomy] -.-> E
    A[Raw FASTQ<br/>96 scats, 2×150 bp<br/>CH16 + CE16 + EU18 pooled] --> B[01 · cutadapt<br/>linked primers at both ends<br/>keeps only CH16 pairs]
    B --> C[02 · DADA2<br/>filter · error model · denoise<br/>merge with trimOverhang]
    C --> D[Chimera removal<br/>+ insert length 60–100 bp]
    D --> E[assignTaxonomy<br/>GenBank chordate 16S<br/>minBoot 80 / 50]
    E --> F[02b · Gulf of California filter<br/>checklist for Actinopteri<br/>remove host and contaminants]
    F --> G[Detection threshold<br/>≥0.1% prey reads and ≥10 reads]
    G --> H[FO · PERMANOVA · PERMDISP<br/>subsampling · rarefaction · Spearman]
    G --> I[04 · Sensitivity<br/>≥1 vs ≥10 reads]
```

## Amplicon and sequencing

| | Sequence (5'→3') |
|---|---|
| Forward primer (CH16_F) | `CGAGAAGACCCTRTGGAGCT` |
| Reverse primer (CH16_R) | `CCTNGGTCGCCCCAAC` |
| Insert without primers | ~77 bp (main peak 76–79 bp; secondary peak 94 bp, *Ophichthus*) |

- Reads are 2 × 150 bp from an **Illumina HiSeq X** (patterned flow cell, binned quality scores) with dual Nextera XT indices.
- Each library pools three amplicons (chordate **CH16**, cephalopod **CE16**, eukaryote **EU18**). The chordate pairs are recovered by their primers.

> [!IMPORTANT]
> The insert is much shorter than the reads, so every read runs through it into the **opposite primer and the Nextera adapter**. Removing only the 5' primers leaves ~100 bp of primer and adapter in each ASV, which inflates ASV counts and distorts taxonomic assignment. Step 1 removes both.

## Reference database

The taxonomy reference (`16S_chor_taxonomy.fasta`) is a global chordate 16S database built from GenBank. It is **not** restricted to the regional checklist. The scripts live in `scripts/00_reference_database/`:

| Step | Script | What it does |
|---|---|---|
| 1 | `01_download_genbank.R` | Downloads mitochondrial records from GenBank with `rentrez`, in batches of 2,000 with retries. The query, record count and date are written to a log. Set your NCBI key in the `ENTREZ_KEY` environment variable. The published database used the query in the script, downloaded in early September 2024. In practice it returned almost every record with "mitochondrion" in the title, including some non-chordates, which the in-silico PCR then filters. A chordate-only query is provided for new builds. |
| 2 | `02_virtual_pcr.R` | Runs an in-silico PCR with the Chord_16S primers (`insect::virtualPCR`, primers kept). |
| 3 | `03_add_taxonomy.R` | Parses accession and species from each title. Removes `UNVERIFIED` records, tRNA entries, unnamed taxa (`sp.`), *Homo sapiens* and *Equus*. Retrieves NCBI taxonomy with `taxize::classification` in resumable batches and writes the 7-rank DADA2 FASTA (`Kingdom;Phylum;Class;Order;Family;Genus;Genus species`). |
| 4 | `04_checklist_coverage.R` | Measures how much of the Gulf of California fish checklist the reference covers, and counts sequences by class. |

```bash
export ENTREZ_KEY=<your NCBI API key>
Rscript scripts/00_reference_database/01_download_genbank.R reference/genbank_chordata_mt.fasta
Rscript scripts/00_reference_database/02_virtual_pcr.R reference/genbank_chordata_mt.fasta reference/chord16S_amplicons.fasta
Rscript scripts/00_reference_database/03_add_taxonomy.R reference/chord16S_amplicons.fasta reference/16S_chor_taxonomy.fasta reference/taxonomy_batches
Rscript scripts/00_reference_database/04_checklist_coverage.R reference/16S_chor_taxonomy.fasta metadata/gulf_of_california_fish_checklist.csv reference_database/
```

**The database used in the study** has 142,689 sequences, 21,686 species and 6,397 genera. By class: Amphibia 39,903, Actinopteri 36,873, Mammalia 28,523, Lepidosauria 16,115, Aves 7,188, Chondrichthyes 3,176, and others (`reference_database/reference_sequences_by_class.csv`).

**Coverage of the Gulf of California fish checklist:**

| Level | In checklist | With a reference sequence | % |
|---|---|---|---|
| Species | 833 | 332 | 39.9 |
| Genus | 432 | 345 | 79.9 |

The per-species list is in `reference_database/checklist_coverage_by_species.csv`. Checklist species without a reference sequence can be misassigned to a congener; this is why the regional filter below truncates such labels to genus.

> [!NOTE]
> Steps 1–3 depend on GenBank and NCBI taxonomy at the time of download (September 2024 for the published database; reference finalised on 22 September 2024). Re-running them today will produce a different database. The exact database used in the study will be archived on Zenodo (*DOI to be added*).

## Gulf of California fish filter

`scripts/02b_gulf_checklist_filter.R` is applied after `assignTaxonomy`:

- **Actinopteri:** a species-level label is kept only if the species is in the Gulf of California fish checklist (`metadata/gulf_of_california_fish_checklist.csv`). Otherwise the label is truncated to genus (`Genus sp`).
- **Checklist source:** the Gulf of California species list of [Shorefishes of the Tropical Eastern Pacific](https://biogeodb.stri.si.edu/sftep/en/pages) (Robertson & Allen 2024, version 3.0, Smithsonian Tropical Research Institute), downloaded in September 2024. It has 833 species and 432 genera and focuses on shore and coastal fishes.
- **Other classes:** not filtered by the checklist.
- **Non-prey taxa are flagged:**
  - host: *Zalophus*;
  - humans and domestic vertebrates: *Homo*, *Bos*, *Gallus*, *Sus*, *Canis*;
  - freshwater: *Oreochromis*.
- **No genus:** ASVs without a genus assignment are excluded from the prey tables.

The script writes an audit table (`results/CH_taxonomy_minBoot80_curated.csv`) with the assigned label, the curated label and the action taken for every ASV.

```bash
Rscript scripts/02b_gulf_checklist_filter.R results/CH_taxonomy_minBoot80.csv \
  metadata/gulf_of_california_fish_checklist.csv results/CH_taxonomy_minBoot80_curated.csv
```

Result at minBoot 80 (2,693 ASVs):

| Action | ASVs |
|---|---|
| Kept at species | 562 |
| Truncated to genus (not in checklist) | 679 |
| Genus only | 322 |
| Excluded: no genus | 977 |
| Non-prey (host, human, domestic) | 152 |
| Non-prey (freshwater) | 1 |

The largest truncations are *Mulloidichthys vanicolensis* (307 ASVs), *Trichiurus lepturus* (86), *Aulopus filamentosus* (54), *Triphoturus mexicanus* (48), *Paranthias colonus* (37) and *Sardinops caeruleus* (35). The same rules are applied inside `03_diet_analysis.R`.

## Quick start

```bash
git clone https://github.com/taniavaldiviac/sealion-diet-metabarcoding.git
cd sealion-diet-metabarcoding

# full run: raw FASTQ folder, reference FASTA, threads
./run_pipeline.sh /path/to/raw_fastqs /path/to/16S_chor_taxonomy.fasta 8
```

To run the steps one at a time:

```bash
# 1. primer and adapter trimming
#    (raw files named CH-<sample>_S<n>_L001_R[12]_001.fastq.gz)
CUTADAPT=cutadapt ./scripts/01_trim_primers_from_raw.sh raw_fastqs/ work/trimmed/ logs/cutadapt/ 8

# 2. DADA2 and taxonomy (writes results/ and logs/)
Rscript scripts/02_dada2.R work/trimmed/ work/filtered/ . reference/16S_chor_taxonomy.fasta 8

# 3. diet analyses (4th argument: folder with earlier outputs to compare, or NA)
Rscript scripts/03_diet_analysis.R . metadata/sample_metadata.csv metadata/gulf_of_california_fish_checklist.csv NA 999

# 4. sensitivity: >= 1 read vs >= 10 reads
Rscript scripts/04_compare_read_thresholds.R analysis/ metadata/sample_metadata.csv
```

## Step-by-step parameters

| Step | Script | Key parameters |
|---|---|---|
| 1 · Trimming | `scripts/01_trim_primers_from_raw.sh` | Linked adapters `^F...rc(R)` on R1 and `^R...rc(F)` on R2, `-e 0.1`, `--discard-untrimmed --pair-filter=any`, `-m 40` |
| 2 · Filtering | `scripts/02_dada2.R` | No truncation, `maxEE = c(2,2)`, `truncQ = 2`, `maxN = 0`, `minLen = 50`, `rm.phix = TRUE` |
| 2 · Error model | `scripts/02_dada2.R` | `learnErrors` with a loess fit forced to be monotonic (for binned quality scores) |
| 2 · Denoise and merge | `scripts/02_dada2.R` | Per-sample `dada`; `mergePairs(minOverlap = 20, maxMismatch = 0, trimOverhang = TRUE)` |
| 2 · Chimeras and length | `scripts/02_dada2.R` | `removeBimeraDenovo(method = "consensus")`; insert 60–100 bp |
| 2 · Taxonomy | `scripts/02_dada2.R` | `assignTaxonomy(tryRC = TRUE, outputBootstraps = TRUE)`; exported at minBoot **80** (main) and 50 |
| 2b · Regional filter | `scripts/02b_gulf_checklist_filter.R` | Gulf of California checklist for Actinopteri; flags non-prey taxa; writes an audit table |
| 3 · Diet analyses | `scripts/03_diet_analysis.R` | Jaccard, 999 permutations, prey present in ≥ 5% of samples (also unfiltered); 999 balanced subsamples; rarefaction to 4,000 prey reads |
| 4 · Sensitivity | `scripts/04_compare_read_thresholds.R` | Presence with ≥ 1 read vs ≥ 10 reads (minBoot 80) |
| extra | `scripts/extra/01b_trim_3prime_from_5prime_trimmed.sh` | Alternative to step 1 for reads already trimmed at 5' (used to cross-check step 1) |

## Curation and detection rules

- 🐟 **Checklist:** species-level assignments within Actinopteri that are not in the Gulf of California fish checklist ([SFTEP](https://biogeodb.stri.si.edu/sftep/en/pages); Robertson & Allen 2024) are reported at genus level (`Genus sp`).
- 🚫 **Non-prey:** removed taxa are the host (*Zalophus*), humans, domestic vertebrates (*Bos*, *Gallus*, *Sus*, *Canis*) and freshwater taxa (*Oreochromis*). ASVs without a genus assignment are excluded from the prey tables.
- 📏 **Detection thresholds**, applied per sample:

  | Threshold | Rule | Role |
  |---|---|---|
  | original | ≥ 1 read | as in the first version of the manuscript |
  | **strict** | ≥ 0.1% of the sample's prey reads **and** ≥ 10 reads | **main result** |

- 📊 **Comparisons:** two one-factor tests.
  - **Interannual:** post-reproductive 2021 (n = 8) vs post-reproductive 2022 (n = 23).
  - **Seasonal:** reproductive 2022 (n = 65) vs post-reproductive 2022 (n = 23).

  PERMDISP uses `betadisper` + `permutest`. Balanced subsampling draws the larger group down to the size of the smaller one in each of 999 iterations.

## Repository structure

```text
sealion-diet-metabarcoding/
├── run_pipeline.sh                 # runs steps 1–4
├── scripts/
│   ├── 00_reference_database/      # GenBank download, in-silico PCR, NCBI taxonomy, checklist coverage
│   ├── 01_trim_primers_from_raw.sh
│   ├── 02_dada2.R
│   ├── 02b_gulf_checklist_filter.R
│   ├── 03_diet_analysis.R
│   ├── 04_compare_read_thresholds.R
│   └── extra/01b_trim_3prime_from_5prime_trimmed.sh
├── metadata/
│   ├── sample_metadata.csv         # year, reproductive season, site
│   ├── gulf_of_california_fish_checklist.csv
│   └── primer_data.csv
├── reference_database/             # checklist coverage and sequences per class
├── results/                        # ASV table, hash key, taxonomy (minBoot 80/50), bootstraps, curated taxonomy
├── analysis/                       # prey matrices, FO, PERMANOVA/PERMDISP, subsampling, Spearman, figures
└── logs/                           # read tracking, parameters, quality and error plots, sessionInfo
```

In `analysis/`, files whose names start with `original_` re-analyse the first (2024) pipeline output, for comparison.

## Software

| Tool | Version |
|---|---|
| cutadapt | 4.4 |
| R | 4.3.3 |
| dada2 | 1.28.0 |
| ShortRead | 1.58.0 |
| vegan | 2.6-10 |
| dplyr · tidyr · readr · tibble | 1.1.4 · 1.3.1 · 2.1.5 · 3.2.1 |
| digest | 0.6.37 |
| rentrez · insect · taxize · Biostrings (reference database) | see `sessionInfo` of the machine that builds it |

Full session information is in `logs/sessionInfo.txt` and `analysis/sessionInfo.txt`.

## Data availability

| Data | Where |
|---|---|
| Raw reads | SRA/ENA: *accession to be added* |
| Reference database | Built with `scripts/00_reference_database/` (see [Reference database](#reference-database)); 142,689 sequences, 21,686 species. Zenodo: *DOI to be added* |

## How to cite

If you use this pipeline, please cite the article (citation to be updated on publication):

```bibtex
@article{mejia_sealion_diet,
  title   = {Fecal eDNA metabarcoding recovers documented temporal dietary patterns
             in California sea lions (Zalophus californianus)},
  author  = {Mej{\'i}a-Guarnizo, Brenda Natasha and others},
  year    = {in review}
}
```

## References

- Deagle B. E., Kirkwood R., Jarman S. N. (2009). Analysis of Australian fur seal diet by pyrosequencing prey DNA in faeces. *Molecular Ecology* 18: 2022–2038.
- Robertson D. R. & Allen G. R. (2024). Shorefishes of the Tropical Eastern Pacific: online information system. Version 3.0. Smithsonian Tropical Research Institute, Balboa, Panamá. https://biogeodb.stri.si.edu/sftep/en/pages (accessed September 2024).

## License

Code is released under the [MIT License](LICENSE).
