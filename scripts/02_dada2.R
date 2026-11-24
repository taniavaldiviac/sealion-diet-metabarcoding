## Step 2: DADA2 ASV inference for the CH16 (Chordata 16S) amplicon.
## Input: FASTQs with both primers and the adapter read-through already removed (step 1).
## Usage: Rscript 02_dada2.R <trimmed_dir> <filtered_dir> <out_dir> <reference_fasta> [threads]

suppressMessages({
  library(dada2); library(ShortRead); library(digest); library(dplyr); library(readr)
  library(tibble); library(tidyr); library(magrittr)
})
args <- commandArgs(trailingOnly = TRUE)
trim_dir <- args[1]; filt_dir <- args[2]; out_dir <- args[3]; ref_fasta <- args[4]
threads <- if (length(args) >= 5) as.integer(args[5]) else 8
set.seed(100)
dir.create(filt_dir, showWarnings = FALSE, recursive = TRUE)
for (d in c("logs", "results")) dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(format(Sys.time(), "%F %T"), "-", ..., "\n")

## Parameters (approved 2026-11-23; minBoot 80 per response letter) -------------------------------------------
PAR <- list(maxEE = c(2, 2), truncQ = 2, maxN = 0, minLen = 50, truncLen = 0,
            minOverlap = 20, maxMismatch = 0, trimOverhang = TRUE,
            chimera_method = "consensus", min_insert = 60, max_insert = 100, minBoot = c(50, 80))  # 80 = main (response letter), 50 = original
writeLines(paste(names(PAR), sapply(PAR, paste, collapse = ","), sep = " = "),
           file.path(out_dir, "logs", "parameters.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "logs", "sessionInfo.txt"))

## Files ----------------------------------------------------------------------
fnFs <- sort(list.files(trim_dir, "_R1_001.fastq.gz$", full.names = TRUE))
fnRs <- sub("_R1_001", "_R2_001", fnFs)
sample.names <- sub("_L001_R1_001.fastq.gz$", "", basename(fnFs))   # e.g. CH-zc78-78-d1-1_S78
filtFs <- file.path(filt_dir, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_dir, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- names(filtRs) <- sample.names
msg(length(fnFs), "samples")

## Quality profiles -----------------------------------------------------------
idx <- seq_len(min(12, length(fnFs)))
ggplot2::ggsave(file.path(out_dir, "logs", "quality_R1.png"), plotQualityProfile(fnFs[idx]), width = 14, height = 9)
ggplot2::ggsave(file.path(out_dir, "logs", "quality_R2.png"), plotQualityProfile(fnRs[idx]), width = 14, height = 9)

## Filter ---------------------------------------------------------------------
msg("filterAndTrim")
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = PAR$truncLen, maxN = PAR$maxN,
                     maxEE = PAR$maxEE, truncQ = PAR$truncQ, minLen = PAR$minLen,
                     rm.phix = TRUE, compress = TRUE, multithread = threads, matchIDs = TRUE)
rownames(out) <- sample.names
keep <- file.exists(filtFs) & file.exists(filtRs) & out[, 2] > 0
filtFs <- filtFs[keep]; filtRs <- filtRs[keep]

## Error model: loess with enforced monotonic decrease (binned quality scores) --
loessErrfun_mono <- function(trans) {
  err <- loessErrfun(trans)
  qs <- as.integer(colnames(err))
  for (i in seq_len(nrow(err))) {
    nm <- rownames(err)[i]
    if (substr(nm, 1, 1) == substr(nm, 3, 3)) next          # skip self-transitions
    err[i, ] <- rev(cummax(rev(err[i, ])))                   # non-increasing with quality
  }
  # recompute self-transitions so each row block sums to 1
  for (nt in c("A", "C", "G", "T")) {
    rows <- grep(paste0("^", nt, "2"), rownames(err))
    self <- paste0(nt, "2", nt)
    err[self, ] <- 1 - colSums(err[setdiff(rownames(err)[rows], self), , drop = FALSE])
  }
  err
}
msg("learnErrors")
errF <- learnErrors(filtFs, errorEstimationFunction = loessErrfun_mono, multithread = threads, randomize = TRUE)
errR <- learnErrors(filtRs, errorEstimationFunction = loessErrfun_mono, multithread = threads, randomize = TRUE)
ggplot2::ggsave(file.path(out_dir, "logs", "errors_R1.png"), plotErrors(errF, nominalQ = TRUE), width = 10, height = 8)
ggplot2::ggsave(file.path(out_dir, "logs", "errors_R2.png"), plotErrors(errR, nominalQ = TRUE), width = 10, height = 8)

## Denoise + merge, sample by sample (low memory) ------------------------------
msg("dada + mergePairs")
mergers <- vector("list", length(filtFs)); names(mergers) <- names(filtFs)
nF <- nR <- setNames(integer(length(filtFs)), names(filtFs))
for (s in names(filtFs)) {
  dF <- derepFastq(filtFs[[s]]); dR <- derepFastq(filtRs[[s]])
  ddF <- dada(dF, err = errF, multithread = threads, verbose = FALSE)
  ddR <- dada(dR, err = errR, multithread = threads, verbose = FALSE)
  nF[s] <- sum(getUniques(ddF)); nR[s] <- sum(getUniques(ddR))
  mergers[[s]] <- mergePairs(ddF, dF, ddR, dR, minOverlap = PAR$minOverlap,
                             maxMismatch = PAR$maxMismatch, trimOverhang = PAR$trimOverhang)
}
seqtab <- makeSequenceTable(mergers)
msg("raw ASVs:", ncol(seqtab))
write_csv(as.data.frame(table(nchar(getSequences(seqtab)))) %>% setNames(c("length", "n_ASV")),
          file.path(out_dir, "logs", "asv_length_distribution_premerge_filters.csv"))

## Chimeras + length filter ----------------------------------------------------
seqtab.nochim <- removeBimeraDenovo(seqtab, method = PAR$chimera_method, multithread = threads)
len <- nchar(colnames(seqtab.nochim))
seqtab.final <- seqtab.nochim[, len >= PAR$min_insert & len <= PAR$max_insert, drop = FALSE]
msg("non-chimeric ASVs:", ncol(seqtab.nochim), "| after length filter:", ncol(seqtab.final))

track <- cbind(out, denoisedF = NA, denoisedR = NA, merged = NA, nonchim = NA, length_filter = NA)
colnames(track)[1:2] <- c("input", "filtered")
s <- names(filtFs)
track[s, "denoisedF"] <- nF; track[s, "denoisedR"] <- nR
track[s, "merged"] <- sapply(mergers, function(m) sum(m$abundance))
track[rownames(seqtab.nochim), "nonchim"] <- rowSums(seqtab.nochim)
track[rownames(seqtab.final), "length_filter"] <- rowSums(seqtab.final)
write_csv(rownames_to_column(as.data.frame(track), "Sample_name"), file.path(out_dir, "logs", "read_tracking.csv"))

## Hashes, ASV table ---------------------------------------------------------
seqs <- colnames(seqtab.final)
hash <- vapply(seqs, digest, "", algo = "sha1", serialize = FALSE)
key <- tibble(Hash = hash, Sequence = seqs, Length = nchar(seqs))
write_csv(key, file.path(out_dir, "results", "CH_hash_key.csv"))
writeLines(paste0(">", key$Hash, "\n", key$Sequence), file.path(out_dir, "results", "CH_hash_key.fasta"))
asv_long <- as.data.frame(seqtab.final) %>% rownames_to_column("Sample_name") %>%
  pivot_longer(-Sample_name, names_to = "Sequence", values_to = "nReads") %>%
  filter(nReads > 0) %>% mutate(Label = "CH", Hash = hash[Sequence]) %>%
  select(Sample_name, Label, Hash, nReads)
write_csv(asv_long, file.path(out_dir, "results", "CH_ASV_table.csv"))
save(seqtab, seqtab.nochim, seqtab.final, track, errF, errR, file = file.path(filt_dir, "dada2_objects.RData"))

## Taxonomy ------------------------------------------------------------------
## One run at the lower threshold; the stricter (main) threshold is derived from the
## stored bootstraps, so both come from the same classification.
msg("assignTaxonomy")
tx <- assignTaxonomy(seqtab.final, ref_fasta, minBoot = min(PAR$minBoot), tryRC = TRUE,
                     outputBootstraps = TRUE, multithread = threads, verbose = TRUE)
levs <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax <- as.data.frame(tx$tax); colnames(tax) <- levs[seq_len(ncol(tax))]
boot <- as.data.frame(tx$boot); colnames(boot) <- levs[seq_len(ncol(boot))]
write_csv(bind_cols(key[, "Hash"], boot), file.path(out_dir, "results", "CH_tax_bootstrap.csv"))
for (mb in PAR$minBoot) {
  t_mb <- tax
  for (j in seq_len(ncol(t_mb))) t_mb[boot[[j]] < mb, j:ncol(t_mb)] <- NA
  write_csv(bind_cols(key, t_mb), file.path(out_dir, "results", paste0("CH_taxonomy_minBoot", mb, ".csv")))
}
msg("done")
