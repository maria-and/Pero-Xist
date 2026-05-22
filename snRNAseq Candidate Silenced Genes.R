# ============================================================
# Final Figure — RNA + H3K27me3 | Xist+ vs Xist-
# FIXED: italic genes + italic Xist + single legend
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)

OUT_DIR <- "~/xci_outputs"
paired_all <- read.csv("~/xci_outputs/paired_all.csv")

# ============================================================
# 1 — DEFINE GENES
# ============================================================

target_genes <- c(
  "Armcx4",
  "Armcx2",
  "Armcx1",
  "Lancl3",
  "Capn6"
)

groups_all <- c(
  "Mus Male",
  "Mus Female",
  "Pero Male",
  "Pero Female"
)

# ============================================================
# 2 — KEEP GENES PRESENT IN BOTH SPECIES
# ============================================================

orth_genes <- paired_all %>%
  filter(display_name %in% target_genes) %>%
  group_by(display_name) %>%
  summarize(n_species = n_distinct(dataset), .groups = "drop") %>%
  filter(n_species == 2) %>%
  pull(display_name)

# ============================================================
# 3 — BUILD DATASET
# ============================================================

gene_df <- paired_all %>%
  filter(display_name %in% orth_genes) %>%
  mutate(group_label = paste(dataset, group)) %>%
  complete(
    display_name,
    group_label = groups_all,
    fill = list(
      rna_pos = 0,
      rna_neg = 0,
      h3k_pos = 0,
      h3k_neg = 0
    )
  ) %>%
  mutate(
    display_name = factor(display_name, levels = target_genes),
    group_label  = factor(group_label, levels = groups_all),
    
    rna_pos = ifelse(is.na(rna_pos), 0, rna_pos),
    rna_neg = ifelse(is.na(rna_neg), 0, rna_neg),
    h3k_pos = ifelse(is.na(h3k_pos), 0, h3k_pos),
    h3k_neg = ifelse(is.na(h3k_neg), 0, h3k_neg),
    
    rna_log2cpm_pos = log2(rna_pos * 1e4 + 1),
    rna_log2cpm_neg = log2(rna_neg * 1e4 + 1),
    h3k_log2cpm_pos = log2(h3k_pos + 1),
    h3k_log2cpm_neg = log2(h3k_neg + 1)
  ) %>%
  select(
    display_name,
    group_label,
    rna_log2cpm_pos,
    rna_log2cpm_neg,
    h3k_log2cpm_pos,
    h3k_log2cpm_neg
  )

# ============================================================
# 4 — LONG FORMAT
# ============================================================

rna_long <- gene_df %>%
  pivot_longer(
    cols = c(rna_log2cpm_pos, rna_log2cpm_neg),
    names_to = "xist_group",
    values_to = "log2cpm"
  ) %>%
  mutate(
    assay = "RNA",
    xist_group = recode(
      xist_group,
      "rna_log2cpm_pos" = "Xist+",
      "rna_log2cpm_neg" = "Xist-"
    )
  )

h3k_long <- gene_df %>%
  pivot_longer(
    cols = c(h3k_log2cpm_pos, h3k_log2cpm_neg),
    names_to = "xist_group",
    values_to = "log2cpm"
  ) %>%
  mutate(
    assay = "H3K27me3",
    xist_group = recode(
      xist_group,
      "h3k_log2cpm_pos" = "Xist+",
      "h3k_log2cpm_neg" = "Xist-"
    )
  )

plot_long <- bind_rows(rna_long, h3k_long) %>%
  mutate(
    xist_group = factor(xist_group, levels = c("Xist+", "Xist-")),
    assay = factor(assay, levels = c("RNA", "H3K27me3"))
  )

# ============================================================
# 5 — COLORS
# ============================================================

rna_colors <- c("Xist+" = "#2166AC", "Xist-" = "#92C5DE")
h3k_colors <- c("Xist+" = "#B2182B", "Xist-" = "#F4A582")

# ============================================================
# 6 — LEGEND LABELS (italic Xist only)
# ============================================================

xist_labels <- c(
  expression(italic(Xist)^"+"),
  expression(italic(Xist)^"-")
)

# ============================================================
# 7 — PANEL FUNCTION
# ============================================================

make_panel <- function(df, group_name, assay_name,
                       color_map, y_lab, tag_letter) {
  
  df_sub <- df %>%
    filter(group_label == group_name, assay == assay_name)
  
  y_max <- max(df_sub$log2cpm, na.rm = TRUE)
  if (!is.finite(y_max)) y_max <- 1
  y_max <- y_max * 1.25
  
  ggplot(df_sub,
         aes(x = display_name, y = log2cpm, fill = xist_group)) +
    
    geom_col(
      position = position_dodge(width = 0.7),
      width = 0.65,
      color = "black",
      linewidth = 0.25
    ) +
    
    geom_text(
      aes(label = sprintf("%.1f", log2cpm)),
      position = position_dodge(width = 0.7),
      vjust = -0.3,
      size = 2.2,
      color = "grey25"
    ) +
    
    scale_fill_manual(
      values = color_map,
      labels = xist_labels
    ) +
    
    # ✔ ROBUST italic gene names
    scale_x_discrete(
      labels = function(x) parse(text = paste0("italic(", x, ")"))
    ) +
    
    coord_cartesian(ylim = c(0, y_max)) +
    
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.05)),
      labels = number_format(accuracy = 0.1)
    ) +
    
    labs(
      title = group_name,
      tag = tag_letter,
      x = NULL,
      y = y_lab,
      fill = NULL
    ) +
    
    theme_classic(base_size = 9) +
    
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.5
      ),
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        size = 8
      ),
      axis.text.y = element_text(size = 7.5),
      legend.position = "none"
    )
}

# ============================================================
# 8 — BUILD PANELS
# ============================================================

groups <- groups_all

rna_panels <- mapply(
  function(g, t) {
    make_panel(plot_long, g, "RNA",
               rna_colors,
               expression("RNA log"[2] * "CPM"), t)
  },
  groups,
  c("A", "B", "C", "D"),
  SIMPLIFY = FALSE
)

h3k_panels <- mapply(
  function(g, t) {
    make_panel(plot_long, g, "H3K27me3",
               h3k_colors,
               expression("H3K27me3 log"[2] * "CPM"), t)
  },
  groups,
  c("E", "F", "G", "H"),
  SIMPLIFY = FALSE
)

# ============================================================
# 9 — SHARED LEGEND
# ============================================================

legend_plot <- make_panel(
  plot_long,
  "Mus Female",
  "RNA",
  rna_colors,
  expression("RNA log"[2] * "CPM"),
  NULL
)

shared_legend <- cowplot::get_legend(
  legend_plot + theme(legend.position = "top")
)

# ============================================================
# 10 — COMBINE FIGURE
# ============================================================

rna_row <- plot_grid(plotlist = rna_panels, ncol = 4)
h3k_row <- plot_grid(plotlist = h3k_panels, ncol = 4)

fig_body <- plot_grid(rna_row, h3k_row, ncol = 1)

final_fig <- plot_grid(
  shared_legend,
  fig_body,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

final_fig

# ============================================================
# 11 — SAVE
# ============================================================

ggsave(
  file.path(OUT_DIR, "final_fixed.png"),
  final_fig,
  width = 14,
  height = 8,
  dpi = 300
)

cat("Saved final_fixed\n")