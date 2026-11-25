## Reference database, step 2: in-silico PCR with the Chord_16S primers (insect::virtualPCR).
## Keeps the primer-binding sites (trimprimers = FALSE), as in the 2024 database.
## Usage: Rscript 02_virtual_pcr.R <genbank_fasta> <out_fasta> [batch_size]
suppressMessages({ library(insect); library(ape) })
args <- commandArgs(trailingOnly = TRUE)
in_fasta <- args[1]; out_fasta <- args[2]
batch_size <- if (length(args) >= 3) as.integer(args[3]) else 1000
CH_F <- "CGAGAAGACCCTRTGGAGCT"
CH_R <- "CCTNGGTCGCCCCAAC"

seqs <- read.dna(in_fasta, format = "fasta")
cat("Input sequences:", length(seqs), "\n")
batches <- split(seq_along(seqs), (seq_along(seqs) - 1) %/% batch_size)
matched <- list()
for (i in seq_along(batches)) {
  res <- tryCatch(virtualPCR(seqs[batches[[i]]], up = CH_F, down = CH_R, trimprimers = FALSE),
                  error = function(e) { warning("Batch ", i, ": ", conditionMessage(e)); NULL })
  if (!is.null(res) && length(res) > 0) matched[[length(matched) + 1]] <- res
  if (i %% 50 == 0) cat("Batch", i, "of", length(batches), "\n")
}
amp <- do.call(c, matched)
cat("Sequences with the amplicon:", length(amp), "\n")
cat(file = out_fasta, paste(paste0(">", names(amp)),
    vapply(dna2char(amp), paste, "", collapse = ""), sep = "\n"), sep = "\n")
