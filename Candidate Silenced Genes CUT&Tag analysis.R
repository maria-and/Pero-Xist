library(dplyr)
library(readr)
library(tidyr)

# ============================================================================
# STEP 1: LOAD DATA
# ============================================================================

cat("=== LOADING DATA ===\n")

mus_h3k27_raw <- read_csv("/Volumes/Extreme Pro/CutandTag_XiXa/mus_rpkm_with_chromosome.csv")
pero_h3k27_raw <- read_csv("/Volumes/Extreme Pro/CutandTag_XiXa/pero_rpkm_with_chromosome.csv")
mus_tpm_raw    <- read_csv("~/Desktop/Docs/musgenome-tpm.csv")
pero_tpm_raw   <- read_csv("~/Desktop/Docs/perogenome-tpm.csv")

# ============================================================================
# STEP 2: PREPARE DATA
# ============================================================================

mus_h3k27 <- mus_h3k27_raw %>%
  select(
    gene_id,
    Mus_H3K27_M1 = Mus_M_1,
    Mus_H3K27_M2 = Mus_M_2,
    Mus_H3K27_F1 = Mus_F_1,
    Mus_H3K27_F2 = Mus_F_2,
    chromosome = chromosome.x
  )

pero_h3k27 <- pero_h3k27_raw %>%
  select(
    gene_id,
    Per_H3K27_M1 = Peromyscus_M_1,
    Per_H3K27_M2 = Peromyscus_M_2,
    Per_H3K27_F1 = Peromyscus_F_1,
    Per_H3K27_F2 = Peromyscus_F_2
  )

mus_tpm <- mus_tpm_raw %>%
  select(
    gene_id = Name,
    Mus_TPM_M1 = M1,
    Mus_TPM_M2 = M2,
    Mus_TPM_F1 = F1,
    Mus_TPM_F2 = F2,
    Mus_TPM_F3 = F3
  )

pero_tpm <- pero_tpm_raw %>%
  select(
    gene_id = Name,
    Per_TPM_M1 = M1,
    Per_TPM_M2 = M2,
    Per_TPM_F1 = F1,
    Per_TPM_F2 = F2,
    Per_TPM_F3 = F3
  )

# ============================================================================
# STEP 3: MERGE ALL DATA
# ============================================================================

combined_clean <- mus_h3k27 %>%
  full_join(pero_h3k27, by = "gene_id") %>%
  full_join(mus_tpm,    by = "gene_id") %>%
  full_join(pero_tpm,   by = "gene_id") %>%
  select(Gene = gene_id, everything()) %>%
  filter(!is.na(Gene))

# ============================================================================
# STEP 4: FILTER X-CHROMOSOME GENES
# ============================================================================

x_genes <- combined_clean %>%
  filter(chromosome %in% c("X", "chrX", "NC_051083.1"))

cat(sprintf("X-linked genes: %d\n", nrow(x_genes)))

# ============================================================================
# SHARED FUNCTIONS
# ============================================================================

perform_ttest <- function(a1, a2, b1, b2, alternative = "greater") {
  g1 <- c(a1, a2)[!is.na(c(a1, a2))]
  g2 <- c(b1, b2)[!is.na(c(b1, b2))]
  if (length(g1) < 2 || length(g2) < 2) return(NA)
  tryCatch(t.test(g1, g2, alternative = alternative)$p.value,
           error = function(e) NA)
}

fisher_p <- function(p1, p2) {
  if (is.na(p1) || is.na(p2)) return(NA)
  pchisq(-2 * (log(p1) + log(p2)), df = 4, lower.tail = FALSE)
}

# ============================================================================
# STEP 5: MALE ANALYSIS
# ============================================================================

male_data <- x_genes %>%
  select(
    Gene,
    Mus_H3K27_M1, Mus_H3K27_M2,
    Per_H3K27_M1, Per_H3K27_M2,
    Mus_TPM_M1, Mus_TPM_M2,
    Per_TPM_M1, Per_TPM_M2
  ) %>%
  filter(!(is.na(Mus_H3K27_M1) & is.na(Per_H3K27_M1))) %>%
  mutate(across(-Gene, ~log2(. + 0.0001))) %>%
  rowwise() %>%
  mutate(
    Mus_H3K27_mean = mean(c(Mus_H3K27_M1, Mus_H3K27_M2), na.rm = TRUE),
    Per_H3K27_mean = mean(c(Per_H3K27_M1, Per_H3K27_M2), na.rm = TRUE),
    Mus_TPM_mean   = mean(c(Mus_TPM_M1, Mus_TPM_M2), na.rm = TRUE),
    Per_TPM_mean   = mean(c(Per_TPM_M1, Per_TPM_M2), na.rm = TRUE),
    
    H3K27_fc = Per_H3K27_mean - Mus_H3K27_mean,
    TPM_fc   = Mus_TPM_mean - Per_TPM_mean,
    combined_score = H3K27_fc + TPM_fc,
    
    H3K27_pval = perform_ttest(Per_H3K27_M1, Per_H3K27_M2,
                               Mus_H3K27_M1, Mus_H3K27_M2),
    TPM_pval   = perform_ttest(Mus_TPM_M1, Mus_TPM_M2,
                               Per_TPM_M1, Per_TPM_M2),
    combined_pval = fisher_p(H3K27_pval, TPM_pval)
  ) %>%
  ungroup()

top_genes_males <- male_data %>%
  filter(!is.na(combined_pval)) %>%
  arrange(combined_pval, desc(combined_score)) %>%
  slice_head(n = 50)

write.csv(male_data, "~/desktop/males_all_xlinked_genes_with_pvalues.csv", row.names = FALSE)
write.csv(top_genes_males, "~/desktop/males_top50_significant_genes.csv", row.names = FALSE)

# ============================================================================
# STEP 6: FEMALE ANALYSIS
# ============================================================================

female_data <- x_genes %>%
  select(
    Gene,
    Mus_H3K27_F1, Mus_H3K27_F2,
    Per_H3K27_F1, Per_H3K27_F2,
    Mus_TPM_F1, Mus_TPM_F2,
    Per_TPM_F1, Per_TPM_F2
  ) %>%
  filter(!(is.na(Mus_H3K27_F1) & is.na(Per_H3K27_F1))) %>%
  mutate(across(-Gene, ~log2(. + 0.1))) %>%
  rowwise() %>%
  mutate(
    Mus_H3K27_mean = mean(c(Mus_H3K27_F1, Mus_H3K27_F2), na.rm = TRUE),
    Per_H3K27_mean = mean(c(Per_H3K27_F1, Per_H3K27_F2), na.rm = TRUE),
    Mus_TPM_mean   = mean(c(Mus_TPM_F1, Mus_TPM_F2), na.rm = TRUE),
    Per_TPM_mean   = mean(c(Per_TPM_F1, Per_TPM_F2), na.rm = TRUE),
    
    H3K27_fc = Per_H3K27_mean - Mus_H3K27_mean,
    TPM_fc   = Mus_TPM_mean - Per_TPM_mean,
    combined_score = H3K27_fc + TPM_fc,
    
    H3K27_pval = perform_ttest(Per_H3K27_F1, Per_H3K27_F2,
                               Mus_H3K27_F1, Mus_H3K27_F2),
    TPM_pval   = perform_ttest(Mus_TPM_F1, Mus_TPM_F2,
                               Per_TPM_F1, Per_TPM_F2),
    combined_pval = fisher_p(H3K27_pval, TPM_pval)
  ) %>%
  ungroup()

top_genes_females <- female_data %>%
  filter(!is.na(combined_pval)) %>%
  arrange(combined_pval, desc(combined_score)) %>%
  slice_head(n = 2000)

write.csv(female_data, "~/desktop/females_all_xlinked_genes_with_pvalues.csv", row.names = FALSE)
write.csv(top_genes_females, "~/desktop/females_all_significant_genes.csv", row.names = FALSE)

# ============================================================================
# STEP 7: MALE vs FEMALE OVERLAP
# ============================================================================

overlap_df <- data.frame(
  Category = c("Both", "Male_only", "Female_only"),
  Count = c(
    length(intersect(top_genes_males$Gene, top_genes_females$Gene)),
    length(setdiff(top_genes_males$Gene, top_genes_females$Gene)),
    length(setdiff(top_genes_females$Gene, top_genes_males$Gene))
  )
)

write.csv(overlap_df, "male_vs_female_gene_overlap.csv", row.names = FALSE)

cat("\nANALYSIS COMPLETE — raw p-values only\n")
