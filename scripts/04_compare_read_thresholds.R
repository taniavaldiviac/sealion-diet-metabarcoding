## Compare prey presence with >=1 read vs >=10 reads (new pipeline, minBoot 80).
## Same PERMANOVA/PERMDISP design as 03_diet_analysis.R (Jaccard, prey in >=5 of 96 samples, 999 perms).
suppressMessages({ library(dplyr); library(readr); library(vegan); library(tibble) })
args <- commandArgs(trailingOnly = TRUE); an <- args[1]; md_file <- args[2]
set.seed(42)
md <- read_csv(md_file, show_col_types = FALSE) %>% rename_with(trimws) %>% rename(Sample_ID = 1) %>%
  filter(grepl("^CH-", Sample_ID)) %>%
  mutate(group = case_when(sampling_year == 2021 ~ "Post2021",
                           reproductive_season == "Posreproductiva" ~ "Post2022", TRUE ~ "Repro2022"))
m <- as.matrix(read.csv(file.path(an, "prey_reads_new_minBoot80.csv"), row.names = 1, check.names = FALSE))
pa <- list(ge1 = (m >= 1) * 1, ge10 = (m >= 10) * 1)
fo <- lapply(pa, function(x) { x <- x[rowSums(x) > 0, ]; colSums(x) / nrow(x) * 100 })
tab <- tibble(Species = colnames(m), FO_ge1 = round(fo$ge1[colnames(m)], 1), FO_ge10 = round(fo$ge10[colnames(m)], 1),
              n_ge1 = colSums(pa$ge1), n_ge10 = colSums(pa$ge10), reads = colSums(m)) %>%
  mutate(FO_ge10 = ifelse(n_ge10 == 0, 0, FO_ge10), diff = FO_ge10 - FO_ge1) %>% arrange(desc(FO_ge1))
write_csv(tab, file.path(an, "FO_ge1_vs_ge10_minBoot80.csv"))
samp <- tibble(Sample_ID = rownames(m), prey_reads = rowSums(m), richness_ge1 = rowSums(pa$ge1), richness_ge10 = rowSums(pa$ge10))
write_csv(samp, file.path(an, "sample_richness_ge1_vs_ge10_minBoot80.csv"))
res <- list()
for (th in names(pa)) for (flt in c("prev5", "none")) {
  x <- pa[[th]]; if (flt == "prev5") x <- x[, colSums(x) >= 5, drop = FALSE]; x <- x[rowSums(x) > 0, , drop = FALSE]
  for (cmp in list(c("Post2021", "Post2022"), c("Repro2022", "Post2022"))) {
    g0 <- md$group[match(rownames(x), md$Sample_ID)]; k <- g0 %in% cmp; g <- factor(g0[k]); xs <- x[k, colSums(x[k, ]) > 0]
    d <- vegdist(xs, "jaccard", binary = TRUE); a <- adonis2(d ~ g, permutations = 999)
    pd <- permutest(betadisper(d, g), permutations = 999)
    res[[length(res) + 1]] <- tibble(threshold = th, filter = flt, comparison = paste(cmp, collapse = " vs "),
      n = paste(table(g)[cmp], collapse = "/"), n_prey = ncol(xs), F = round(a$F[1], 2), R2 = round(a$R2[1], 3),
      p = a$`Pr(>F)`[1], disp_F = round(pd$tab$F[1], 2), disp_p = pd$tab$`Pr(>F)`[1])
  }
}
write_csv(bind_rows(res), file.path(an, "PERMANOVA_ge1_vs_ge10_minBoot80.csv"))
print(as.data.frame(bind_rows(res)))
cat("prey taxa ge1:", sum(tab$n_ge1 > 0), " ge10:", sum(tab$n_ge10 > 0), "\n")
cat("empty samples ge1:", sum(samp$richness_ge1 == 0), " ge10:", sum(samp$richness_ge10 == 0), "\n")
print(as.data.frame(samp %>% filter(grepl("zc06-|zc12-|zc41-", Sample_ID))))
print(as.data.frame(tab %>% filter(n_ge10 == 0)))
print(as.data.frame(head(tab, 20)))
print(summary(samp$richness_ge1 - samp$richness_ge10))
s <- cor.test(samp$richness_ge1, samp$prey_reads, method = "spearman", exact = FALSE); s2 <- cor.test(samp$richness_ge10, samp$prey_reads, method = "spearman", exact = FALSE)
cat("spearman ge1", s$estimate, s$p.value, " ge10", s2$estimate, s2$p.value, "\n")
