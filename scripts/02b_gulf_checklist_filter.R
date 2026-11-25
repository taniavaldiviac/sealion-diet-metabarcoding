## Step 2b: Gulf of California fish filter and removal of non-prey taxa.
## Bony-fish (Actinopteri) species-level labels are kept only if the species is in the regional
## checklist (Robertson et al. 2024); otherwise the assignment is truncated to genus ("Genus sp").
## Other classes are not filtered by the checklist. Host, human, domestic and freshwater taxa are
## flagged as non-prey. The same rules are applied inside 03_diet_analysis.R; this script writes
## the audit table so every change can be inspected.
## Usage: Rscript 02b_gulf_checklist_filter.R <taxonomy_csv> <checklist_csv> <out_csv>
suppressMessages({ library(dplyr); library(readr); library(stringr) })
args <- commandArgs(trailingOnly = TRUE)
tax_file <- args[1]; checklist_file <- args[2]; out_file <- args[3]
NON_PREY_GENERA   <- c("Zalophus", "Homo", "Bos", "Gallus", "Sus", "Canis")  # host, human, domestic
FRESHWATER_GENERA <- c("Oreochromis")                                       # tilapia
LEVS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

gulf <- read_csv(checklist_file, show_col_types = FALSE) %>% rename(species = 1) %>%
  mutate(species = str_squish(species)) %>% filter(!is.na(species), species != "") %>% pull(species)
tax <- read_csv(tax_file, show_col_types = FALSE) %>% mutate(across(any_of(LEVS), ~ na_if(.x, "NA")))

out <- tax %>% mutate(
  Species_assigned = Species,
  genus_any = coalesce(Genus, word(Species, 1)),
  action = case_when(
    genus_any %in% NON_PREY_GENERA ~ "non-prey: host/human/domestic",
    genus_any %in% FRESHWATER_GENERA ~ "non-prey: freshwater",
    is.na(Genus) & is.na(Species) ~ "excluded: no genus assignment",
    !is.na(Class) & Class == "Actinopteri" & !is.na(Species) & !(Species %in% gulf) ~ "truncated to genus: not in Gulf checklist",
    is.na(Species) ~ "genus only (Genus sp)",
    TRUE ~ "kept at species"),
  Species = case_when(
    action == "truncated to genus: not in Gulf checklist" ~ paste(Genus, "sp"),
    action == "genus only (Genus sp)" ~ paste(Genus, "sp"),
    TRUE ~ Species),
  prey = !str_starts(action, "non-prey") & action != "excluded: no genus assignment") %>%
  select(-genus_any)
write_csv(out, out_file)
print(as.data.frame(count(out, action, name = "ASVs")))
cat("\nSpecies truncated to genus (not in the Gulf checklist):\n")
print(as.data.frame(out %>% filter(action == "truncated to genus: not in Gulf checklist") %>%
  count(Species_assigned, Species, name = "ASVs") %>% arrange(desc(ASVs))), row.names = FALSE)
