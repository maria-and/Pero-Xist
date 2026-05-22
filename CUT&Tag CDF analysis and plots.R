###############################################
# Load libraries
###############################################
library(tidyverse)
library(patchwork)
library(scales)

###############################################
# Load your files
###############################################
mus_rpkm_df  <- read_csv("/Volumes/Extreme Pro/CutandTag_XiXa/mus_rpkm_with_chromosome.csv")
pero_rpkm_df <- read_csv("/Volumes/Extreme Pro/CutandTag_XiXa/pero_rpkm.csv")

mus_female_cols <- c("Mus_F_1", "Mus_F_2")
mus_male_cols   <- c("Mus_M_1", "Mus_M_2")
pero_female_cols <- c("Pero_F_1", "Pero_F_2")
pero_male_cols   <- c("Pero_M_1", "Pero_M_2")

###############################################
# Prepare ECDF data
###############################################
prepare_ecdf_data <- function(df, male_cols, female_cols, gene_col){
  df %>%
    filter(!is.na(Chromosome),
           Chromosome != "",
           Chromosome != "MT",
           Chromosome != "M") %>%
    mutate(chr_type = ifelse(Chromosome == "X", "X", "Autosome")) %>%
    filter(
      rowSums(select(., all_of(male_cols)) > 0) > 0 &
        rowSums(select(., all_of(female_cols)) > 0) > 0
    ) %>%
    mutate(
      mean_rpkm_male   = rowMeans(select(., all_of(male_cols)), na.rm = TRUE),
      mean_rpkm_female = rowMeans(select(., all_of(female_cols)), na.rm = TRUE),
      log2_male   = log2(mean_rpkm_male + 0.001),
      log2_female = log2(mean_rpkm_female + 0.001)
    ) %>%
    select(all_of(gene_col), Chromosome, chr_type, log2_male, log2_female)
}

mus_prepped <- prepare_ecdf_data(mus_rpkm_df, mus_male_cols, mus_female_cols, "gene_id")
pero_prepped <- prepare_ecdf_data(pero_rpkm_df, pero_male_cols, pero_female_cols, "Name")

split_data <- function(df){
  list(
    x = df %>% filter(Chromosome == "X"),
    autosome = df %>% filter(Chromosome != "X")
  )
}

mus_split  <- split_data(mus_prepped)
pero_split <- split_data(pero_prepped)

###############################################
# Check sample sizes
###############################################
cat("\n=== Sample Sizes ===\n")
cat("Mus X genes:       ", nrow(mus_split$x), "\n")
cat("Mus Autosome genes:", nrow(mus_split$autosome), "\n")
cat("Pero X genes:      ", nrow(pero_split$x), "\n")
cat("Pero Autosome genes:", nrow(pero_split$autosome), "\n\n")

###############################################
# Get X chromosome sample sizes for matching
###############################################
x_sample_size <- min(nrow(mus_split$x), nrow(pero_split$x))
cat("X chromosome sample size:", x_sample_size, "\n\n")

# Use SMALLER sample size for cross-species comparisons
cross_species_subsample <- round(x_sample_size )  # 50% of X chr size
cat("Cross-species autosome subsample:", cross_species_subsample, "\n\n")

###############################################
# KS Test function with optional subsampling
###############################################
adaptive_ks_test <- function(data1, data2, subsample_n = NULL) {
  set.seed(42)  # For reproducibility
  
  data1 <- na.omit(data1)
  data2 <- na.omit(data2)
  
  if (!is.null(subsample_n)) {
    n1 <- min(length(data1), subsample_n)
    n2 <- min(length(data2), subsample_n)
    data1 <- sample(data1, n1, replace = FALSE)
    data2 <- sample(data2, n2, replace = FALSE)
  }
  
  # Run KS test
  ks_result <- ks.test(data1, data2)
  return(ks_result$p.value)
}

###############################################
# Run all tests - use smaller samples for cross-species
###############################################
###############################################
# Set seed for reproducible results
###############################################
set.seed(123)  # Use any number you want

###############################################
# Run all tests - will now be reproducible
###############################################

# 1. Pl female vs male
pl_fem_male <- list(
  x = adaptive_ks_test(pero_split$x$log2_female, pero_split$x$log2_male),
  auto = adaptive_ks_test(pero_split$autosome$log2_female, 
                          pero_split$autosome$log2_male, 
                          subsample_n = x_sample_size)
)

# 2. Pl female vs Mm female
pl_mm_female <- list(
  x = adaptive_ks_test(pero_split$x$log2_female, mus_split$x$log2_female),
  auto = adaptive_ks_test(pero_split$autosome$log2_female, 
                          mus_split$autosome$log2_female,
                          subsample_n = x_sample_size)
)

# 3. Pl female vs Mm male
pl_f_mm_m <- list(
  x = adaptive_ks_test(pero_split$x$log2_female, mus_split$x$log2_male),
  auto = adaptive_ks_test(pero_split$autosome$log2_female, 
                          mus_split$autosome$log2_male,
                          subsample_n = x_sample_size)
)

# 4. Mm female vs male
mm_fem_male <- list(
  x = adaptive_ks_test(mus_split$x$log2_female, mus_split$x$log2_male),
  auto = adaptive_ks_test(mus_split$autosome$log2_female, 
                          mus_split$autosome$log2_male, 
                          subsample_n = x_sample_size)
)

# 5. Pl male vs Mm male
pl_mm_male <- list(
  x = adaptive_ks_test(pero_split$x$log2_male, mus_split$x$log2_male),
  auto = adaptive_ks_test(pero_split$autosome$log2_male, 
                          mus_split$autosome$log2_male,
                          subsample_n = x_sample_size)
)

###############################################
# Print results
###############################################
results_df <- data.frame(
  Comparison = c("Pl female vs male", "Pl female vs Mm female", 
                 "Pl female vs Mm male", "Mm female vs male",
                 "Pl male vs Mm male"),
  X_Chromosome = c(pl_fem_male$x, pl_mm_female$x, 
                   pl_f_mm_m$x, mm_fem_male$x,
                   pl_mm_male$x),
  Autosomes = c(pl_fem_male$auto, pl_mm_female$auto,
                pl_f_mm_m$auto, mm_fem_male$auto,
                pl_mm_male$auto)
)

cat("=== Results (KS Test, adaptive subsampling) ===\n")
print(results_df, digits = 3)
###############################################
# Build combined ECDF data for all groups
###############################################
build_combined_ecdf <- function() {
  
  # X chromosome data
  x_data <- data.frame(
    bin_mid = fine_bins,
    Pl_F = ecdf(pero_split$x$log2_female)(fine_bins),
    Pl_M = ecdf(pero_split$x$log2_male)(fine_bins),
    Mm_F = ecdf(mus_split$x$log2_female)(fine_bins),
    Mm_M = ecdf(mus_split$x$log2_male)(fine_bins)
  )
  
  # Autosome data
  auto_data <- data.frame(
    bin_mid = fine_bins,
    Pl_F = ecdf(pero_split$autosome$log2_female)(fine_bins),
    Pl_M = ecdf(pero_split$autosome$log2_male)(fine_bins),
    Mm_F = ecdf(mus_split$autosome$log2_female)(fine_bins),
    Mm_M = ecdf(mus_split$autosome$log2_male)(fine_bins)
  )
  
  list(x = x_data, auto = auto_data)
}

ecdf_data <- build_combined_ecdf()

###############################################
# Pivot to long format for ggplot
###############################################
x_long <- ecdf_data$x %>%
  pivot_longer(cols = -bin_mid, names_to = "Group", values_to = "CumulativeFreq")

auto_long <- ecdf_data$auto %>%
  pivot_longer(cols = -bin_mid, names_to = "Group", values_to = "CumulativeFreq")

###############################################
# Define colors and labels
###############################################
color_palette <- c(
  "Pl_F" = "red",    # Red
  "Pl_M" = "darkblue",    # Blue
  "Mm_F" = "pink",    # Orange
  "Mm_M" = "lightblue"     # Purple
)

group_labels <- c(
  "Pl_F" = "P. leucopus Female",
  "Pl_M" = "P. leucopus Male",
  "Mm_F" = "M. musculus Female",
  "Mm_M" = "M. musculus Male"
)

###############################################
# Create X chromosome plot
###############################################
p_x <- ggplot(x_long, aes(x = bin_mid, y = CumulativeFreq, color = Group)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  scale_color_manual(
    values = color_palette,
    labels = group_labels,
    name = NULL
  ) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    x = "Log2(RPKM)",
    y = "Cumulative Frequency",
    title = "X Chromosome"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2))

###############################################
# Create Autosome plot
###############################################
p_auto <- ggplot(auto_long, aes(x = bin_mid, y = CumulativeFreq, color = Group)) +
  geom_line(linewidth = 1.2, alpha = 0.8) +
  scale_color_manual(
    values = color_palette,
    labels = group_labels,
    name = NULL
  ) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  labs(
    x = "Log2(RPKM)",
    y = "Cumulative Frequency",
    title = "Autosomes"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    legend.position = "bottom",
    legend.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(nrow = 2))

###############################################
# Combine into final 2-panel plot
###############################################
final_plot <- p_x | p_auto

final_plot

# Save
ggsave("expression_comparison_2panel.pdf", final_plot, width = 14, height = 6)
ggsave("expression_comparison_2panel.png", final_plot, width = 14, height = 6, dpi = 300)

