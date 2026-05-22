# USED FOR FIG5B UPDATED ON 03092026

# ============================
# Load libraries
# ============================
library(tidyverse)
library(pheatmap)
library(cowplot)
library(grid)
library(RColorBrewer)
library(ggplotify)

# ============================
# Load data
# ============================
df <- read.csv("~/Desktop/Docs/males_all_xlinked_genes_with_pvalues.csv")

# ============================
# p-value → stars
# ============================
p_to_stars <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ ""
  )
}

# ============================
# Combined significance ONLY (Fisher)
# ============================
df <- df %>%
  mutate(
    combined_pval = pchisq(
      -2 * log(TPM_pval * H3K27_pval),
      df = 4,
      lower.tail = FALSE
    ),
    combined_sig = p_to_stars(combined_pval)
  )

# ============================
# Color palette
# ============================
heat_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(200)

# ============================
# PANEL A — RNA ranked by RNA (TPM)
# ============================
df_rna <- df %>%
  arrange(desc(TPM_fc)) %>%
  slice_head(n = 50)

rna_mat <- df_rna %>%
  select(Mus_TPM_mean, Per_TPM_mean) %>%
  as.matrix()

rownames(rna_mat) <- paste0(df_rna$Gene, "  ", df_rna$combined_sig)

rna_plot <- pheatmap(
  rna_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = heat_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "Gene Expression",
  border_color = NA,
  legend = TRUE,
  silent = TRUE
)

# ============================
# PANEL B — H3 ranked by H3 signal
# ============================
df_h3 <- df %>%
  arrange(desc(H3K27_fc)) %>%
  slice_head(n = 50)

h3_mat <- df_h3 %>%
  select(Mus_H3K27_mean, Per_H3K27_mean) %>%
  as.matrix()

rownames(h3_mat) <- paste0(df_h3$Gene, "  ", df_h3$combined_sig)

h3_plot <- pheatmap(
  h3_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  scale = "none",
  color = heat_colors,
  fontsize_row = 8,
  fontsize_col = 8,
  main = "H3K27me3 Signal",
  border_color = NA,
  legend = TRUE,
  silent = TRUE
)

# ============================
# Combine plots
# ============================
combined1 <- plot_grid(
  as.grob(rna_plot$gtable),
  as.grob(h3_plot$gtable),
  ncol = 2,
  align = "h",
  label_size = 12
)

print(combined1)

