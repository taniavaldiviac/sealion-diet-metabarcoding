## Step 3: taxonomic curation, frequency of occurrence (FO), PERMANOVA, PERMDISP,
## balanced subsampling, rarefaction check and richness-depth correlation.
## The same analysis is applied to every dataset (new pipeline at minBoot 80/50 and the
## original 2024 pipeline re-thresholded at 50/80) so the results are directly comparable.
## Usage: Rscript 03_diet_analysis.R <out_dir> <metadata_csv> <gulf_fish_checklist_csv> [original_results_dir|NA] [n_iter]
##
## Design (manuscript section 2.4 and response letter):
##  - Prey reads exclude the host (Zalophus californianus) and contaminants
##    (human, Bos, Gallus, Sus, freshwater taxa listed in FRESHWATER_GENERA).
##  - Detection thresholds: "original" = >=1 read; "strict" = >=0.1% of the sample's prey
##    reads AND >=10 reads.
##  - PERMANOVA on presence/absence, Jaccard, 999 permutations, two one-factor comparisons:
##      interannual: post-reproductive 2021 vs post-reproductive 2022
##      seasonal:    reproductive 2022 vs post-reproductive 2022
##    with and without removing prey present in <5% of all samples.
##  - PERMDISP (betadisper) for each comparison.
##  - Balanced subsampling: n_iter iterations, larger group subsampled to the smaller one.
##  - Rarefaction to 4,000 prey reads per sample as a check.
##  - Spearman correlation between per-sample prey richness and prey read depth.

suppressMessages({ library(dplyr); library(tidyr); library(readr); library(vegan); library(tibble); library(ggplot2) })
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[1]; md_file <- args[2]; checklist_file <- args[3]; orig_dir <- args[4]
N_ITER <- if (length(args) >= 5) as.integer(args[5]) else 999
res_dir <- file.path(out_dir, "results"); an_dir <- file.path(out_dir, "analysis")
dir.create(an_dir, showWarnings = FALSE)
set.seed(42)
NPERM <- 999; RAREFY_DEPTH <- 4000; MIN_PREV <- 0.05
LEVS <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
THRESHOLDS <- list(original = c(rel = 0, reads = 1), strict = c(rel = 0.001, reads = 10))
NON_PREY_GENERA <- c("Zalophus", "Homo", "Bos", "Gallus", "Sus", "Canis")
FRESHWATER_GENERA <- c("Oreochromis")   # from analysis/taxa_to_review.csv (tilapia, freshwater)

## Metadata --------------------------------------------------------------------
md <- read_csv(md_file, show_col_types = FALSE) %>% rename_with(trimws) %>% rename(Sample_ID = 1) %>%
  filter(grepl("^CH-", Sample_ID)) %>%
  mutate(group = case_when(sampling_year == 2021 & reproductive_season == "Posreproductiva" ~ "Post2021",
                           sampling_year == 2022 & reproductive_season == "Posreproductiva" ~ "Post2022",
                           sampling_year == 2022 & reproductive_season == "Reproductiva" ~ "Repro2022"))
COMPARISONS <- list(interannual = c("Post2021", "Post2022"), seasonal = c("Repro2022", "Post2022"))
gulf <- read_csv(checklist_file, show_col_types = FALSE) %>% rename(Species = 1) %>% pull(Species) %>% trimws()

## Curation (rules of the original report, Report_lobo_html.qmd) ------------------
curate <- function(tax) {
  tax %>% mutate(across(all_of(LEVS), ~ na_if(.x, "NA")),
    Species = ifelse(!is.na(Class) & Class == "Actinopteri" & !is.na(Species) & !(Species %in% gulf), NA, Species),
    Species = ifelse(is.na(Species) & !is.na(Genus), paste(Genus, "sp"), Species),
    Species = ifelse(Species %in% c("Eukaryotic synthetic"), NA, Species),
    prey = !is.na(Species) &
      !(coalesce(Genus, sub(" .*", "", Species)) %in% c(NON_PREY_GENERA, FRESHWATER_GENERA)))
}

resolution <- function(tax, asv, label) {
  reads <- asv %>% group_by(Hash) %>% summarise(reads = sum(nReads), .groups = "drop")
  t <- tax %>% mutate(across(all_of(LEVS), ~ na_if(.x, "NA"))) %>% inner_join(reads, by = "Hash")
  bind_rows(lapply(c("Family", "Genus", "Species"), function(l) tibble(dataset = label, level = l,
    pct_ASVs = round(100 * mean(!is.na(t[[l]])), 1),
    pct_reads = round(100 * sum(t$reads[!is.na(t[[l]])]) / sum(t$reads), 1))))
}

prey_matrix <- function(tc, asv) {
  m <- asv %>% inner_join(filter(tc, prey) %>% select(Hash, Species), by = "Hash") %>%
    mutate(Sample_ID = sub("_S[0-9]+$", "", Sample_name)) %>%
    group_by(Sample_ID, Species) %>% summarise(reads = sum(nReads), .groups = "drop") %>%
    pivot_wider(names_from = Species, values_from = reads, values_fill = 0) %>%
    column_to_rownames("Sample_ID") %>% as.matrix()
  m[intersect(md$Sample_ID, rownames(m)), , drop = FALSE]
}

presence <- function(m, th) {
  rel <- m / pmax(rowSums(m), 1)
  pa <- (m >= th["reads"] & rel >= th["rel"]) * 1
  pa[rowSums(pa) > 0, colSums(pa) > 0, drop = FALSE]
}

fo_table <- function(pa, m) {
  g <- md$group[match(rownames(pa), md$Sample_ID)]
  base <- tibble(Species = colnames(pa), n_total = colSums(pa), FO_total = round(100 * colMeans(pa), 1),
                 reads = colSums(m[rownames(pa), colnames(pa), drop = FALSE]))
  for (gr in sort(unique(g))) base[[paste0("FO_", gr)]] <- round(100 * colMeans(pa[g == gr, , drop = FALSE]), 1)
  arrange(base, desc(FO_total))
}

one_test <- function(pa_sub, g) {
  pa_sub <- pa_sub[, colSums(pa_sub) > 0, drop = FALSE]
  d <- vegdist(pa_sub, method = "jaccard", binary = TRUE)
  a <- adonis2(d ~ g, permutations = NPERM)
  pd <- permutest(betadisper(d, g), permutations = NPERM)
  c(F = a$F[1], R2 = a$R2[1], p = a$`Pr(>F)`[1], disp_F = pd$tab$F[1], disp_p = pd$tab$`Pr(>F)`[1])
}

multivariate <- function(pa, label, th_name, subsample = FALSE) {
  out <- list(); sub_out <- list()
  for (flt in c("prev5pct", "unfiltered")) {
    x <- if (flt == "prev5pct") pa[, colMeans(pa) >= MIN_PREV, drop = FALSE] else pa
    x <- x[rowSums(x) > 0, , drop = FALSE]
    for (cmp in names(COMPARISONS)) {
      g_all <- md$group[match(rownames(x), md$Sample_ID)]
      xs <- x[g_all %in% COMPARISONS[[cmp]], , drop = FALSE]
      g <- factor(g_all[g_all %in% COMPARISONS[[cmp]]])
      r <- one_test(xs, g)
      out[[length(out) + 1]] <- tibble(dataset = label, threshold = th_name, prevalence_filter = flt,
        comparison = cmp, n = paste(table(g), collapse = "/"), n_prey = sum(colSums(xs) > 0),
        PERMANOVA_F = r["F"], PERMANOVA_R2 = r["R2"], PERMANOVA_p = r["p"],
        PERMDISP_F = r["disp_F"], PERMDISP_p = r["disp_p"])
      if (subsample && flt == "prev5pct") {
        nmin <- min(table(g)); big <- names(which.max(table(g)))
        it <- t(replicate(N_ITER, {
          keep <- c(which(g != big), sample(which(g == big), nmin))
          one_test(xs[keep, , drop = FALSE], droplevels(g[keep]))
        }))
        sub_out[[length(sub_out) + 1]] <- as_tibble(it) %>%
          mutate(dataset = label, threshold = th_name, comparison = cmp, iter = row_number(), .before = 1)
      }
    }
  }
  list(main = bind_rows(out), sub = bind_rows(sub_out))
}

## Datasets -------------------------------------------------------------------
asv_new <- read_csv(file.path(res_dir, "CH_ASV_table.csv"), show_col_types = FALSE)
datasets <- list()
for (mb in c(80, 50)) datasets[[paste0("new_minBoot", mb)]] <- list(asv = asv_new,
  tax = read_csv(file.path(res_dir, paste0("CH_taxonomy_minBoot", mb, ".csv")), show_col_types = FALSE))
## Optional: original 2024 pipeline outputs (1_CH_ASV_table.csv, 1_CH_tax_bootstrap.csv) for comparison
if (!is.na(orig_dir) && dir.exists(orig_dir)) {
asv_orig <- read_csv(file.path(orig_dir, "1_CH_ASV_table.csv"), show_col_types = FALSE)
tb_orig <- read_csv(file.path(orig_dir, "1_CH_tax_bootstrap.csv"), show_col_types = FALSE)
for (mb in c(50, 80)) {
  t <- tb_orig %>% select(Hash, Sequence, all_of(LEVS))
  for (j in seq_along(LEVS)) { b <- tb_orig[[paste0("boot.", LEVS[j])]]; t[!is.na(b) & b < mb, LEVS[j:length(LEVS)]] <- NA }
  datasets[[paste0("original_minBoot", mb)]] <- list(asv = asv_orig, tax = t)
}
}
MAIN <- "new_minBoot80"

## Run ------------------------------------------------------------------------
summ <- res <- fo_all <- mv <- sub <- rar <- spe <- list(); labels <- list()
for (nm in names(datasets)) {
  ds <- datasets[[nm]]; tc <- curate(ds$tax)
  labels[[nm]] <- tc %>% select(Hash, Species, prey)
  res[[nm]] <- resolution(ds$tax, ds$asv, nm)
  m <- prey_matrix(tc, ds$asv)
  write_csv(rownames_to_column(as.data.frame(m), "Sample_ID"), file.path(an_dir, paste0("prey_reads_", nm, ".csv")))
  summ[[nm]] <- tibble(dataset = nm, ASVs = n_distinct(ds$asv$Hash), total_reads = sum(ds$asv$nReads),
    prey_ASVs = sum(tc$prey & tc$Hash %in% ds$asv$Hash), prey_reads = sum(m), samples_with_prey = sum(rowSums(m) > 0),
    prey_taxa_any_read = ncol(m), prey_taxa_strict = ncol(presence(m, THRESHOLDS$strict)))
  for (th in names(THRESHOLDS)) {
    pa <- presence(m, THRESHOLDS[[th]])
    fo_all[[paste(nm, th)]] <- fo_table(pa, m) %>% mutate(dataset = nm, threshold = th, .before = 1)
    r <- multivariate(pa, nm, th, subsample = nm == MAIN)
    mv[[paste(nm, th)]] <- r$main; sub[[paste(nm, th)]] <- r$sub
    depth <- rowSums(m[rownames(pa), , drop = FALSE])
    ct <- cor.test(rowSums(pa), depth, method = "spearman", exact = FALSE)
    spe[[paste(nm, th)]] <- tibble(dataset = nm, threshold = th, n = nrow(pa), rho = ct$estimate, p = ct$p.value)
  }
  ## Rarefaction check (4,000 prey reads; presence = >=1 read after rarefying)
  mm <- m[rowSums(m) >= RAREFY_DEPTH, , drop = FALSE]
  rr <- rrarefy(mm, RAREFY_DEPTH); pa_r <- (rr > 0) * 1
  pa_r <- pa_r[, colSums(pa_r) > 0, drop = FALSE]
  rar[[nm]] <- multivariate(pa_r, nm, "rarefied4000")$main
  fo_all[[paste(nm, "rarefied4000")]] <- fo_table(pa_r, mm) %>% mutate(dataset = nm, threshold = "rarefied4000", .before = 1)
  ct <- cor.test(rowSums(pa_r), rowSums(mm), method = "spearman", exact = FALSE)
  spe[[paste(nm, "rar")]] <- tibble(dataset = nm, threshold = "rarefied4000", n = nrow(pa_r), rho = ct$estimate, p = ct$p.value)
}

write_csv(bind_rows(summ), file.path(an_dir, "summary_by_dataset.csv"))
write_csv(bind_rows(res), file.path(an_dir, "taxonomic_resolution.csv"))
fo <- bind_rows(fo_all); write_csv(fo, file.path(an_dir, "FO_all.csv"))
write_csv(bind_rows(c(mv, rar)), file.path(an_dir, "PERMANOVA_PERMDISP.csv"))
s <- bind_rows(sub); write_csv(s, file.path(an_dir, "balanced_subsampling_iterations.csv"))
s %>% group_by(dataset, threshold, comparison) %>%
  summarise(across(c(F, R2, p, disp_p), list(median = median, q025 = ~ quantile(.x, .025), q975 = ~ quantile(.x, .975))),
            prop_p_lt_0.05 = mean(p < 0.05), .groups = "drop") %>%
  write_csv(file.path(an_dir, "balanced_subsampling_summary.csv"))
write_csv(bind_rows(spe), file.path(an_dir, "richness_depth_spearman.csv"))

## Prey label changes between minBoot 50 and 80 (new pipeline)
reads_h <- asv_new %>% group_by(Hash) %>% summarise(reads = sum(nReads), n_samples = n(), .groups = "drop")
labels$new_minBoot50 %>% rename(label50 = Species, prey50 = prey) %>%
  inner_join(labels$new_minBoot80 %>% rename(label80 = Species, prey80 = prey), by = "Hash") %>%
  inner_join(reads_h, by = "Hash") %>% filter(prey50 | prey80) %>%
  group_by(label50, label80) %>% summarise(ASVs = n(), reads = sum(reads), .groups = "drop") %>%
  mutate(change = case_when(is.na(label80) ~ "lost (no genus at 80)",
                            label50 == label80 ~ "same",
                            grepl(" sp$", label80) & !grepl(" sp$", label50) ~ "species -> genus",
                            TRUE ~ "different")) %>%
  arrange(change, desc(reads)) %>% write_csv(file.path(an_dir, "prey_label_changes_minBoot50_vs_80.csv"))

## Taxa to review as possible contaminants / non-marine (non-fish classes, genera not in checklist)
tc <- curate(datasets[[MAIN]]$tax) %>% filter(prey) %>% inner_join(reads_h, by = "Hash")
tc %>% filter(Class != "Actinopteri" | !(sub(" .*", "", Species) %in% sub(" .*", "", gulf))) %>%
  group_by(Class, Order, Family, Genus, Species) %>% summarise(reads = sum(reads), .groups = "drop") %>%
  arrange(desc(reads)) %>% write_csv(file.path(an_dir, "taxa_to_review.csv"))

## FO > 5% (Figure 8) and NMDS for the main result
fo8 <- fo %>% filter(dataset == MAIN, threshold == "strict", FO_total > 5)
write_csv(fo8, file.path(an_dir, "Fig8_FO_gt5pct_main.csv"))
p <- ggplot(fo8, aes(reorder(Species, FO_total), FO_total)) + geom_col(fill = "grey40") + coord_flip() +
  labs(x = NULL, y = "Frequency of occurrence (%)") + theme_minimal(base_size = 11)
ggsave(file.path(an_dir, "Fig8_FO_main.png"), p, width = 7, height = 0.25 * nrow(fo8) + 1.5, dpi = 300)
m <- as.matrix(read.csv(file.path(an_dir, paste0("prey_reads_", MAIN, ".csv")), row.names = 1, check.names = FALSE))
pa <- presence(m, THRESHOLDS$strict); pa <- pa[, colMeans(pa) >= MIN_PREV, drop = FALSE]; pa <- pa[rowSums(pa) > 0, ]
nm <- metaMDS(vegdist(pa, "jaccard", binary = TRUE), k = 2, trymax = 100, trace = 0)
sc <- as.data.frame(scores(nm)); sc$group <- md$group[match(rownames(sc), md$Sample_ID)]
p <- ggplot(sc, aes(NMDS1, NMDS2, colour = group)) + geom_point(size = 2) + stat_ellipse() +
  labs(subtitle = sprintf("Jaccard, stress = %.3f", nm$stress)) + theme_minimal(base_size = 11)
ggsave(file.path(an_dir, "NMDS_main.png"), p, width = 6, height = 5, dpi = 300)
writeLines(capture.output(sessionInfo()), file.path(an_dir, "sessionInfo.txt"))
cat("done\n")
