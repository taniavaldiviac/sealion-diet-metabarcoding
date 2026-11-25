## Reference database, step 3: clean names, retrieve NCBI taxonomy and write a DADA2 reference.
## Header format for assignTaxonomy: >Kingdom;Phylum;Class;Order;Family;Genus;Genus species
## Usage: Rscript 03_add_taxonomy.R <amplicon_fasta> <out_fasta> [work_dir]
## An NCBI API key is read from ENTREZ_KEY (recommended for taxize).
suppressMessages({ library(Biostrings); library(dplyr); library(stringr); library(taxize); library(purrr); library(tidyr) })
args <- commandArgs(trailingOnly = TRUE)
in_fasta <- args[1]; out_fasta <- args[2]
work <- if (length(args) >= 3) args[3] else "taxonomy_batches"
dir.create(work, showWarnings = FALSE, recursive = TRUE)

## 1. Accession and species name from the GenBank title
x <- readDNAStringSet(in_fasta)
df <- tibble(title = names(x), sequence = as.character(x)) %>%
  mutate(accession = str_extract(title, "^\\S+"),
         raw_name = str_remove(title, "^\\S+\\s+"),
         species = str_extract(raw_name, "^[A-Z][a-z]+ [a-z]+"))

## 2. Remove unusable or unwanted records
df <- df %>%
  filter(!str_detect(raw_name, "UNVERIFIED"),
         !str_detect(raw_name, "Phe-|Val-tRNA"),
         !is.na(species),
         !str_detect(species, " sp$"),
         !species %in% c("Homo sapiens"),          # human removed from the reference
         !str_detect(species, "^Equus "),          # as in the 2024 database
         species != "Indirana semipalmata",
         !str_detect(sequence, "^N+$"))
cat("Records after cleaning:", nrow(df), "| species:", n_distinct(df$species), "\n")

## 3. NCBI taxonomy in batches (cached, so the step can be resumed)
sp <- unique(df$species)
batches <- split(sp, ceiling(seq_along(sp) / 500))
for (i in seq_along(batches)) {
  f <- file.path(work, sprintf("batch_%04d.rds", i))
  if (file.exists(f)) next
  res <- tryCatch(classification(batches[[i]], db = "ncbi", rows = 1, messages = FALSE),
                  error = function(e) { message("Batch ", i, ": ", conditionMessage(e)); NULL })
  if (!is.null(res)) saveRDS(res, f)
}
ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
to_row <- function(r, s) {
  if (!is.data.frame(r) || !all(c("rank", "name") %in% names(r)))
    return(tibble(species = s, !!!setNames(rep(list(NA_character_), 6), ranks)))
  r %>% filter(rank %in% ranks) %>% distinct(rank, .keep_all = TRUE) %>% select(rank, name) %>%
    pivot_wider(names_from = rank, values_from = name) %>% mutate(species = s)
}
tax <- list.files(work, "^batch_.*\\.rds$", full.names = TRUE) %>%
  map_dfr(function(f) imap_dfr(readRDS(f), to_row))
for (r in ranks) if (!r %in% names(tax)) tax[[r]] <- NA_character_
# NCBI uses "Metazoa" at the kingdom rank for animals
write.csv(tax, file.path(work, "taxonomy_table.csv"), row.names = FALSE)

## 4. Write the DADA2 reference
ref <- df %>% inner_join(tax, by = "species") %>%
  mutate(header = paste(kingdom, phylum, class, order, family, genus, species, sep = ";"))
writeLines(paste0(">", ref$header, "\n", ref$sequence), out_fasta)
cat("Reference written:", nrow(ref), "sequences,", n_distinct(ref$species), "species,",
    n_distinct(ref$genus), "genera\n")
