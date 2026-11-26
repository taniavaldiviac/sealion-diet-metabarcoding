## Reference database, step 4: coverage of the Gulf of California fish checklist.
## Checklist source: Shorefishes of the Tropical Eastern Pacific online information system
## (Robertson D. R. & Allen G. R., Smithsonian Tropical Research Institute),
## https://biogeodb.stri.si.edu/sftep/en/pages — species list for the Gulf of California,
## downloaded in September 2024 (metadata/gulf_of_california_fish_checklist.csv).
## Reports which checklist species and genera have at least one reference sequence.
## Usage: Rscript 04_checklist_coverage.R <reference_fasta> <checklist_csv> <out_dir>
suppressMessages({ library(dplyr); library(readr); library(stringr) })
args <- commandArgs(trailingOnly = TRUE)
ref_fasta <- args[1]; checklist <- args[2]; out_dir <- args[3]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

hdr <- sub("^>", "", grep("^>", readLines(ref_fasta), value = TRUE))
ref <- tibble(header = hdr) %>%
  mutate(class = str_split_i(header, ";", 3), genus = str_split_i(header, ";", 6),
         species = str_split_i(header, ";", 7))
ck <- read_csv(checklist, show_col_types = FALSE) %>% rename(species = 1) %>%
  mutate(species = str_squish(species)) %>% filter(!is.na(species), species != "") %>%
  mutate(genus = word(species, 1)) %>% distinct(species, .keep_all = TRUE)

ck <- ck %>% mutate(species_in_reference = species %in% ref$species,
                    genus_in_reference = genus %in% ref$genus)
cov <- tibble(
  level = c("species", "genus"),
  checklist = c(nrow(ck), n_distinct(ck$genus)),
  in_reference = c(sum(ck$species_in_reference), n_distinct(ck$genus[ck$genus_in_reference]))) %>%
  mutate(pct = round(100 * in_reference / checklist, 1))
print(as.data.frame(cov))
write_csv(cov, file.path(out_dir, "checklist_coverage_summary.csv"))
write_csv(ck %>% select(species, genus, species_in_reference, genus_in_reference),
          file.path(out_dir, "checklist_coverage_by_species.csv"))
write_csv(ref %>% count(class, name = "sequences") %>% arrange(desc(sequences)),
          file.path(out_dir, "reference_sequences_by_class.csv"))
cat("Reference:", length(hdr), "sequences,", n_distinct(ref$species), "species,", n_distinct(ref$genus), "genera\n")
