## Reference database, step 1: download chordate mitochondrial sequences from GenBank.
## Usage: Rscript 01_download_genbank.R <out_fasta> ["<entrez query>"]
## An NCBI API key is read from the ENTREZ_KEY environment variable (never hard-code it).
## Default = the query as written in the 2024 script (chordata/rdb/16S_NCBI.R). Note that it
## says 12S[Title] and has an unbalanced parenthesis; mitochondrial genomes enter through
## "mitochondrion[Title]". TO CONFIRM which query and date produced the published database.
DEFAULT_QUERY <- 'Chordata[Organism] AND (12S[Title]) AND (partial[Title] OR complete[Title]) OR mitochondrion[Title])'

suppressMessages(library(rentrez))
args <- commandArgs(trailingOnly = TRUE)
out_fasta <- args[1]
query <- if (length(args) >= 2) args[2] else DEFAULT_QUERY
if (nzchar(Sys.getenv("ENTREZ_KEY"))) set_entrez_key(Sys.getenv("ENTREZ_KEY"))

search <- entrez_search(db = "nucleotide", term = query, use_history = TRUE)
cat("Query:", query, "\nRecords:", search$count, "\nDownloaded on:", format(Sys.Date()), "\n")
chunk <- 2000
con <- file(out_fasta, open = "wt")
for (start in seq(0, search$count - 1, by = chunk)) {
  seqs <- NULL
  for (attempt in 1:5) {                       # retry on transient NCBI errors
    seqs <- tryCatch(entrez_fetch(db = "nucleotide", web_history = search$web_history,
                                  rettype = "fasta", retmax = chunk, retstart = start),
                     error = function(e) NULL)
    if (!is.null(seqs)) break
    Sys.sleep(5 * attempt)
  }
  if (is.null(seqs)) stop("Failed to download records starting at ", start)
  writeLines(seqs, con)
  cat("Records", start + 1, "-", min(start + chunk, search$count), "\n")
  Sys.sleep(0.4)
}
close(con)
writeLines(c(paste("query:", query), paste("records:", search$count), paste("date:", Sys.Date())),
           paste0(out_fasta, ".download_log.txt"))
