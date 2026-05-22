# ============================================================
# XCI Analysis — Peromyscus leucopus + Mus musculus
# Full pipeline: DEG CSVs + GTFs + H3K27me3 (sex-split)
#
# NORMALIZATION STRATEGY (cross-species snPairedTag):
#   H3K: CPM per cell → pseudo-bulk mean per group → log2FC
#        + z-score across X-linked genes within species
#        (removes species-level peak density / sequencing depth differences)
#   RNA: log-normalized counts (Seurat NormalizeData scale.factor=1e4)
#        → pseudo-bulk mean per group → log2FC
#        + z-score across X-linked genes within species
#   For within-species plots: use raw CPM / log-norm means (absolute)
#   For cross-species plots:  use z-scored values
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(rtracklayer)
library(Matrix)
library(patchwork)
library(ggrepel)
library(scales)
library(ggVennDiagram)
library(tidyr)
library(GenomicRanges)

# ============================================================
# CONFIGURATION
# ============================================================

path_pl     <- "/Volumes/Extreme Pro/pairedtag_pero/ShaSun_H3K27me3_pl.rawSeurat.RDS"
path_mm     <- "/Volumes/Extreme Pro/pairedtag_pero/ShaSun_H3K27me3_mm.rawSeurat.RDS"
gtf_path_pl <- "~/Downloads/Peromyscus_numbered_fixed.gtf"
gtf_path_mm <- "~/Downloads/Mus_musculus.GRCm39.115.gtf.gz"

deg_mus_m_f  <- "~/xci_outputs/deg_mus_male_xlinked.csv"
deg_mus_f_f  <- "~/xci_outputs/deg_mus_female_xlinked.csv"
deg_pero_m_f <- "~/xci_outputs/deg_pero_male_xlinked.csv"
deg_pero_f_f <- "~/xci_outputs/deg_pero_female_xlinked.csv"
orthologs_f  <- "~/xci_outputs/orthologs_x_name_matched.csv"

bc_pl_m_pos <- "~/barcodes_pl_male_xist_pos.txt"
bc_pl_m_neg <- "~/barcodes_pl_male_xist_neg.txt"
bc_pl_f_pos <- "~/barcodes_pl_female_xist_pos.txt"
bc_pl_f_neg <- "~/barcodes_pl_female_xist_neg.txt"

bc_mm_m_pos <- "~/barcodes_mm_male_xist_pos.txt"
bc_mm_m_neg <- "~/barcodes_mm_male_xist_neg.txt"
bc_mm_f_pos <- "~/barcodes_mm_female_xist_pos.txt"
bc_mm_f_neg <- "~/barcodes_mm_female_xist_neg.txt"

# NOTE: double-check these — sample names must match sobj@meta.data$sample exactly.
# Run: table(readRDS(path_pl)@meta.data$sample) to verify before running.
PERO_MALE_SAMPLES   <- c("PL_M_1", "PL_F_2")
PERO_FEMALE_SAMPLES <- c("PL_F_1", "PL_M_2")

MUS_MALE_SAMPLES   <- c("6_1_M_R1", "7_1_M_R2")
MUS_FEMALE_SAMPLES <- c("8_1_F_R1", "9_1_F_R2")

OUT_DIR <- "~/xci_outputs/"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

BASE <- 14

pub_theme <- theme_classic(base_size = BASE) +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = BASE,   color = "black"),
    axis.line        = element_line(linewidth = 0.5,              color = "black"),
    axis.ticks       = element_line(linewidth = 0.5,              color = "black"),
    axis.text        = element_text(size = BASE - 1,              color = "black"),
    axis.title       = element_text(size = BASE,                  color = "black"),
    legend.text      = element_text(size = BASE - 1,              color = "black"),
    legend.title     = element_text(size = BASE - 1,              color = "black"),
    legend.key.size  = unit(0.5, "cm"),
    legend.position  = "none",
    panel.spacing    = unit(0.4, "cm"),
    plot.margin      = margin(4, 6, 4, 6, "pt")
  )

cat_colors <- c(
  "Silenced + H3K27me3 gain" = "#B5294E",
  "Silenced only"            = "#E8A0A0",
  "H3K27me3 gain only"       = "#2F6DB5",
  "Other"                    = "grey82"
)

fix_cols <- function(df) {
  colnames(df) <- iconv(colnames(df), to = "ASCII//TRANSLIT")
  df
}

fmt_p <- function(p) {
  if (is.na(p))   return("n.s.")
  if (p < 0.0001) return("****")
  if (p < 0.001)  return("***")
  if (p < 0.01)   return("**")
  if (p < 0.05)   return("*")
  return("n.s.")
}


# ============================================================
# BLOCK 1 — Load DEG CSVs + orthologs
# ============================================================
cat("Loading DEG CSVs...\n")

deg_mus_m  <- fix_cols(read.csv(deg_mus_m_f))
deg_mus_f  <- fix_cols(read.csv(deg_mus_f_f))
deg_pero_m <- fix_cols(read.csv(deg_pero_m_f))
deg_pero_f <- fix_cols(read.csv(deg_pero_f_f))
orthologs  <- read.csv(orthologs_f)

cat("Column names (Mus Male DEG):\n")
print(colnames(deg_mus_m))

rank_silenced <- function(df, group_label) {
  df %>%
    filter(avg_log2FC < 0) %>%
    arrange(avg_log2FC) %>%
    mutate(rank  = row_number(),
           group = group_label)
}

ranked_mus_m  <- rank_silenced(deg_mus_m,  "Mus Male")
ranked_mus_f  <- rank_silenced(deg_mus_f,  "Mus Female")
ranked_pero_m <- rank_silenced(deg_pero_m, "Pero Male")
ranked_pero_f <- rank_silenced(deg_pero_f, "Pero Female")

cat(sprintf(
  "Silenced genes — Mus Male: %d  Mus Female: %d  Pero Male: %d  Pero Female: %d\n",
  nrow(ranked_mus_m), nrow(ranked_mus_f),
  nrow(ranked_pero_m), nrow(ranked_pero_f)
))

all_ranked <- bind_rows(ranked_mus_m, ranked_mus_f,
                        ranked_pero_m, ranked_pero_f)
write.csv(all_ranked,
          file.path(OUT_DIR, "ranked_silenced_all_groups.csv"),
          row.names = FALSE)
cat("Saved ranked_silenced_all_groups.csv\n")


# ============================================================
# BLOCK 2 — Load GTFs + build TSS windows
# ============================================================
cat("\nLoading GTFs...\n")

load_gtf <- function(gtf_path, x_chr) {
  gtf_raw <- as.data.frame(import(gtf_path))
  gtf     <- gtf_raw[gtf_raw$type == "gene", ]
  rm(gtf_raw); gc()
  
  chrom_col <- intersect(c("seqnames","chr","chrom","Chr"), colnames(gtf))[1]
  name_col  <- intersect(c("gene_name","Name","gene"),      colnames(gtf))[1]
  
  if (is.na(chrom_col))
    stop(sprintf("No chromosome column found in GTF: %s", gtf_path))
  if (is.na(name_col))
    warning(sprintf("No gene-name column in GTF %s — using gene_id", gtf_path))
  
  df <- data.frame(
    gene_id    = gtf$gene_id,
    gene_name  = if (!is.na(name_col)) gtf[[name_col]] else gtf$gene_id,
    chromosome = as.character(gtf[[chrom_col]]),
    start      = gtf$start,
    end        = gtf$end,
    stringsAsFactors = FALSE
  )
  df$display_name <- ifelse(
    !is.na(df$gene_name) & df$gene_name != "" & df$gene_name != df$gene_id,
    df$gene_name, df$gene_id
  )
  df$is_x <- df$chromosome == x_chr
  cat(sprintf("  GTF loaded: %d genes, %d on chrX (%s)\n",
              nrow(df), sum(df$is_x), x_chr))
  df
}

gtf_mm <- load_gtf(gtf_path_mm, "X")
gtf_pl <- load_gtf(gtf_path_pl, "X")

tss_windows <- function(gtf_df, window = 2000) {
  gtf_df %>%
    filter(is_x) %>%
    mutate(tss       = start,
           win_start = tss - window,
           win_end   = tss + window)
}

tss_mm <- tss_windows(gtf_mm)
tss_mm$chromosome[tss_mm$chromosome == "X"] <- "chrX"

gtf_pl$chromosome[gtf_pl$chromosome == "X"] <- "NC-051083.1"
gtf_pl$is_x <- gtf_pl$chromosome == "NC-051083.1"
tss_pl <- tss_windows(gtf_pl)

cat(sprintf("TSS windows — Mus: %d  Pero: %d\n", nrow(tss_mm), nrow(tss_pl)))


# ============================================================
# BLOCK 2B — Generate sex-split Xist barcode files
# ============================================================
generate_sex_barcodes <- function(seurat_path, species,
                                  male_samples, female_samples,
                                  pos_out_m, neg_out_m,
                                  pos_out_f, neg_out_f,
                                  xist_gene = NULL,          # <-- ADD THIS
                                  qc_rna   = 750,
                                  qc_peaks = 750,
                                  force    = FALSE) {
  
  all_exist <- all(file.exists(pos_out_m, neg_out_m, pos_out_f, neg_out_f))
  if (all_exist && !force) {
    cat(sprintf("[%s] Barcode files exist — skipping (set force=TRUE to redo)\n",
                species))
    return(invisible(NULL))
  }
  
  # ── 1. Load object ──────────────────────────────────────────
  cat(sprintf("[%s] Reading RDS...\n", species))
  sobj <- readRDS(seurat_path)
  
  # ── 2. Pull ONLY what we need, then free the object ─────────
  meta         <- sobj@meta.data                         # data.frame, tiny
  barcodes_all <- rownames(meta)
  
  DefaultAssay(sobj) <- "RNA"
  rna_mat_all <- GetAssayData(sobj, layer = "counts")   # ALL genes
  
  # Full library size per cell — this is what corrects the depth confound
  col_totals_full <- Matrix::colSums(rna_mat_all)
  col_totals_full[col_totals_full == 0] <- 1
  
  # ── Extract Xist counts ───────────────────────────────────────
  if (is.null(xist_gene)) {
    # fallback: case-insensitive grep
    hits <- grep("^xist$", rownames(rna_mat_all), ignore.case = TRUE, value = TRUE)
    if (length(hits) == 0)
      stop(sprintf("[%s] Xist gene not found. Rows containing 'xist': %s",
                   species,
                   paste(grep("xist", rownames(rna_mat_all),
                              ignore.case = TRUE, value = TRUE), collapse = ", ")))
    xist_gene <- hits[1]
  }
  
  if (!xist_gene %in% rownames(rna_mat_all))
    stop(sprintf("[%s] Specified xist_gene '%s' not found in RNA matrix.",
                 species, xist_gene))
  
  xist_counts <- as.numeric(rna_mat_all[xist_gene, ])
  names(xist_counts) <- colnames(rna_mat_all)
  cat(sprintf("[%s] Xist gene used: '%s'\n", species, xist_gene))
  
  # Subset to X-linked genes only AFTER computing full library size
  n_by_name <- sum(tss_df$gene_name %in% rownames(rna_mat_all))
  n_by_id   <- sum(tss_df$gene_id   %in% rownames(rna_mat_all))
  cat(sprintf("  RNA match — gene_name: %d  |  gene_id: %d\n", n_by_name, n_by_id))
  rna_key_col <- if (n_by_name >= n_by_id) "gene_name" else "gene_id"
  cat(sprintf("  Using '%s' to match RNA rows\n", rna_key_col))
  
  x_keys    <- unique(tss_df[[rna_key_col]][tss_df[[rna_key_col]] %in% rownames(rna_mat_all)])
  rna_mat_x <- rna_mat_all[x_keys, , drop = FALSE]
  rm(rna_mat_all); gc()
  
  # Normalize X-linked counts by FULL library size — depth confound corrected
  rna_cpm <- sweep(rna_mat_x, 2, col_totals_full / 1e6, "/")
  
  
  # ── CRITICAL: free the large object before any further work ──
  rm(sobj); gc(); gc()
  cat(sprintf("[%s] Seurat object freed.\n", species))
  
  # ── 3. QC filter (operates on metadata vectors only) ─────────
  required_qc_cols <- c("nCount_RNA", "nCount_peaks", "sample")
  missing_cols     <- setdiff(required_qc_cols, colnames(meta))
  if (length(missing_cols) > 0)
    stop(sprintf("[%s] Missing metadata columns: %s",
                 species, paste(missing_cols, collapse = ", ")))
  
  keep <- meta$nCount_RNA   >= qc_rna  &
    meta$nCount_peaks >= qc_peaks &
    meta$sample       != "multiplet"
  
  meta         <- meta[keep, ]
  xist_counts  <- xist_counts[keep]
  barcodes_all <- barcodes_all[keep]
  
  cat(sprintf("[%s] Cells after QC: %d\n", species, sum(keep)))
  cat(sprintf("[%s] Sample table:\n", species))
  print(table(meta$sample))
  
  # ── 4. Sex assignment ─────────────────────────────────────────
  sex_derived <- dplyr::case_when(
    meta$sample %in% male_samples   ~ "male",
    meta$sample %in% female_samples ~ "female",
    TRUE                             ~ NA_character_
  )
  
  n_na <- sum(is.na(sex_derived))
  if (n_na == nrow(meta))
    stop(sprintf(
      paste0("[%s] Zero cells matched.\n",
             "  male_samples:   %s\n",
             "  female_samples: %s\n",
             "  Actual sample values (first 8): %s"),
      species,
      paste(male_samples,   collapse = ", "),
      paste(female_samples, collapse = ", "),
      paste(unique(meta$sample)[seq_len(min(8, length(unique(meta$sample))))],
            collapse = ", ")
    ))
  if (n_na > 0)
    warning(sprintf("[%s] %d cells unassigned — verify sample name lists",
                    species, n_na))
  
  cat(sprintf("[%s] Sex assignment:\n", species))
  print(table(sex_derived, useNA = "ifany"))
  
  # ── 5. Xist positivity ───────────────────────────────────────
  xist_pos <- xist_counts > 0
  
  # Depth diagnostic — now xist_pos actually exists
  cat(sprintf("[%s] RNA depth median — Xist+: %.0f   Xist-: %.0f\n",
              species,
              median(col_totals_full[names(col_totals_full) %in% barcodes_all[ xist_pos]]),
              median(col_totals_full[names(col_totals_full) %in% barcodes_all[!xist_pos]])))
  
  # ── 6. Write barcode files ────────────────────────────────────
  is_male   <- !is.na(sex_derived) & sex_derived == "male"
  is_female <- !is.na(sex_derived) & sex_derived == "female"
  
  writeLines(barcodes_all[ is_male   &  xist_pos], pos_out_m)
  writeLines(barcodes_all[ is_male   & !xist_pos], neg_out_m)
  writeLines(barcodes_all[ is_female &  xist_pos], pos_out_f)
  writeLines(barcodes_all[ is_female & !xist_pos], neg_out_f)
  
  cat(sprintf(
    "[%s] Saved:\n  Male   Xist+: %d → %s\n  Male   Xist-: %d → %s\n  Female Xist+: %d → %s\n  Female Xist-: %d → %s\n",
    species,
    sum(is_male   &  xist_pos), pos_out_m,
    sum(is_male   & !xist_pos), neg_out_m,
    sum(is_female &  xist_pos), pos_out_f,
    sum(is_female & !xist_pos), neg_out_f
  ))
  
  rm(meta, xist_counts); gc()
}

# Pero — known Xist gene
generate_sex_barcodes(
  path_pl, "Pero",
  male_samples   = PERO_MALE_SAMPLES,
  female_samples = PERO_FEMALE_SAMPLES,
  pos_out_m = bc_pl_m_pos, neg_out_m = bc_pl_m_neg,
  pos_out_f = bc_pl_f_pos, neg_out_f = bc_pl_f_neg,
  xist_gene = "LOC114708611",   # <-- explicit
  force = TRUE
)

# Mus — fill in after running the grep above
generate_sex_barcodes(
  path_mus, "Mus",
  male_samples   = MUS_MALE_SAMPLES,
  female_samples = MUS_FEMALE_SAMPLES,
  pos_out_m = bc_mus_m_pos, neg_out_m = bc_mus_m_neg,
  pos_out_f = bc_mus_f_pos, neg_out_f = bc_mus_f_neg,
  xist_gene = "snmus_XXXXXXX",  # <-- paste result from grep
  force = TRUE
)
# ============================================================
# BLOCK 3 — Extract H3K27me3 AND RNA signal at TSS windows
#
# NORMALIZATION:
#   H3K: CPM per cell (library-size corrected within cell)
#        → pseudo-bulk mean per group (Xist+ vs Xist-)
#        → log2FC = log2((pos_mean + 1) / (neg_mean + 1))
#        → z_score_h3k: z-score of h3k_log2fc across all X genes
#          within this species+sex combination
#          (makes Pero and Mus values directly comparable)
#
#   RNA: log-normalized per cell (log1p(counts / total * 1e4))
#        → pseudo-bulk mean per group
#        → log2FC
#        → z_score_rna: z-score across X genes within species+sex
#
#   Both raw means (h3k_pos, h3k_neg, rna_pos, rna_neg) are retained
#   for absolute within-species plots.
# ============================================================

extract_paired_xist_groups <- function(seurat_path,
                                       tss_df,
                                       xist_pos_barcodes,
                                       xist_neg_barcodes,
                                       species_label,
                                       sex_label,
                                       qc_rna   = 750,
                                       qc_peaks = 750,
                                       peak_sep = "-") {
  
  cat(sprintf("\n===== [%s %s] =====\n", species_label, sex_label))
  
  # 1. Load + QC
  sobj <- readRDS(seurat_path)
  keep_mask <- sobj$nCount_RNA   >= qc_rna  &
    sobj$nCount_peaks >= qc_peaks &
    sobj$sample       != "multiplet"
  sobj <- sobj[, keep_mask]
  gc()
  cat(sprintf("  After QC: %d cells\n", ncol(sobj)))
  
  # 2. Auto-detect peak assay
  pk_assay <- grep("peak|atac|h3k", Assays(sobj),
                   ignore.case = TRUE, value = TRUE)[1]
  if (is.na(pk_assay))
    stop(sprintf("[%s %s] No peak assay found. Available: %s",
                 species_label, sex_label,
                 paste(Assays(sobj), collapse = ", ")))
  cat(sprintf("  Peak assay detected: '%s'\n", pk_assay))
  
  # 3. Barcode intersection
  all_cells <- colnames(sobj)
  pos_cells <- intersect(all_cells, xist_pos_barcodes)
  neg_cells <- intersect(all_cells, xist_neg_barcodes)
  cat(sprintf("  Xist+ = %d  |  Xist- = %d\n",
              length(pos_cells), length(neg_cells)))
  
  if (length(pos_cells) < 10 || length(neg_cells) < 10)
    stop(sprintf("[%s %s] Too few cells: Xist+ %d, Xist- %d",
                 species_label, sex_label,
                 length(pos_cells), length(neg_cells)))
  
  sobj <- sobj[, union(pos_cells, neg_cells)]
  gc()
  
  # 4. H3K peak matrix — CPM per cell
  DefaultAssay(sobj) <- pk_assay
  peak_mat   <- GetAssayData(sobj, layer = "counts")
  peak_names <- rownames(peak_mat)
  
  if (peak_sep == ":") {
    pk_chr   <- sub(":([0-9]+)-([0-9]+)$",   "", peak_names)
    pk_start <- as.integer(sub(".*:([0-9]+)-[0-9]+$",   "\\1", peak_names))
    pk_end   <- as.integer(sub(".*:([0-9]+)-([0-9]+)$", "\\2", peak_names))
  } else {
    pk_chr   <- sub("-([0-9]+)-([0-9]+)$",   "", peak_names)
    pk_start <- as.integer(sub(".*-([0-9]+)-[0-9]+$",   "\\1", peak_names))
    pk_end   <- as.integer(sub(".*-([0-9]+)$",           "\\1", peak_names))
  }
  
  cat(sprintf("  Peak chromosomes (unique, first 12): %s\n",
              paste(head(sort(unique(pk_chr)), 12), collapse = ", ")))
  cat(sprintf("  TSS  chromosomes (unique):           %s\n",
              paste(sort(unique(tss_df$chromosome)), collapse = ", ")))
  
  chr_overlap <- intersect(unique(pk_chr), unique(tss_df$chromosome))
  cat(sprintf("  Matching chromosomes:                %s\n",
              if (length(chr_overlap) > 0)
                paste(chr_overlap, collapse = ", ")
              else "*** NONE ***"))
  
  if (length(chr_overlap) == 0)
    stop(sprintf("[%s %s] Zero chromosome matches.", species_label, sex_label))
  
  peak_gr <- GRanges(seqnames = pk_chr,
                     ranges   = IRanges(pk_start, pk_end))
  tss_gr  <- GRanges(seqnames = tss_df$chromosome,
                     ranges   = IRanges(tss_df$win_start, tss_df$win_end))
  
  hits <- findOverlaps(tss_gr, peak_gr)
  cat(sprintf("  Peak-TSS overlaps: %d\n", length(hits)))
  if (length(hits) == 0)
    stop(sprintf("[%s %s] findOverlaps returned 0 hits.", species_label, sex_label))
  
  needed_peaks <- sort(unique(subjectHits(hits)))
  peak_mat     <- peak_mat[needed_peaks, , drop = FALSE]
  
  hits_df <- data.frame(gene_idx = queryHits(hits),
                        peak_idx = subjectHits(hits))
  hits_df <- hits_df[hits_df$peak_idx %in% needed_peaks, ]
  hits_df$peak_idx <- match(hits_df$peak_idx, needed_peaks)
  
  col_totals_pk <- Matrix::colSums(peak_mat)
  col_totals_pk[col_totals_pk == 0] <- 1
  pk_cpm <- sweep(peak_mat, 2, col_totals_pk / 1e6, "/")
  rm(peak_mat); gc()
  
  # 5. RNA matrix — log-norm per cell (auto-detect gene_name vs gene_id)
  DefaultAssay(sobj) <- "RNA"
  rna_mat <- GetAssayData(sobj, layer = "counts")
  
  n_by_name <- sum(tss_df$gene_name %in% rownames(rna_mat))
  n_by_id   <- sum(tss_df$gene_id   %in% rownames(rna_mat))
  cat(sprintf("  RNA gene match — by gene_name: %d  |  by gene_id: %d\n",
              n_by_name, n_by_id))
  
  rna_key_col <- if (n_by_name >= n_by_id) "gene_name" else "gene_id"
  cat(sprintf("  Using '%s' to match RNA rows\n", rna_key_col))
  
  x_keys    <- unique(tss_df[[rna_key_col]][tss_df[[rna_key_col]] %in% rownames(rna_mat)])
  rna_mat_x <- rna_mat[x_keys, , drop = FALSE]
  
  col_totals_rna <- Matrix::colSums(rna_mat)
  col_totals_rna[col_totals_rna == 0] <- 1
  rna_lognorm <- log1p(sweep(rna_mat_x, 2, col_totals_rna / 1e4, "/"))
  rm(rna_mat); gc()
  
  # 6. Group indices
  pk_cells    <- colnames(pk_cpm)
  pos_idx     <- which(pk_cells %in% pos_cells)
  neg_idx     <- which(pk_cells %in% neg_cells)
  
  rna_cells   <- colnames(rna_lognorm)
  rna_pos_idx <- which(rna_cells %in% pos_cells)
  rna_neg_idx <- which(rna_cells %in% neg_cells)
  
  # 7. Gene-level aggregation — H3K
  n_genes       <- nrow(tss_df)
  gene_h3k_pos  <- rep(NA_real_, n_genes)
  gene_h3k_neg  <- rep(NA_real_, n_genes)
  gene_h3k_pval <- rep(NA_real_, n_genes)
  hits_by_gene  <- split(hits_df$peak_idx, hits_df$gene_idx)
  
  for (g in seq_len(n_genes)) {
    pk <- hits_by_gene[[as.character(g)]]
    if (is.null(pk)) next
    sub_mat          <- pk_cpm[pk, , drop = FALSE]
    pos_vals         <- Matrix::colSums(sub_mat[, pos_idx, drop = FALSE])
    neg_vals         <- Matrix::colSums(sub_mat[, neg_idx, drop = FALSE])
    gene_h3k_pos[g]  <- mean(pos_vals, na.rm = TRUE)
    gene_h3k_neg[g]  <- mean(neg_vals, na.rm = TRUE)
    if (length(pos_vals) > 3 && length(neg_vals) > 3)
      gene_h3k_pval[g] <- wilcox.test(pos_vals, neg_vals,
                                      exact = FALSE)$p.value
  }
  
  # 8. Gene-level aggregation — RNA (uses rna_key_col)
  gene_rna_pos  <- rep(NA_real_, n_genes)
  gene_rna_neg  <- rep(NA_real_, n_genes)
  gene_rna_pval <- rep(NA_real_, n_genes)
  
  for (g in seq_len(n_genes)) {
    key <- tss_df[[rna_key_col]][g]
    if (!key %in% rownames(rna_lognorm)) next
    pos_vals         <- as.numeric(rna_lognorm[key, rna_pos_idx])
    neg_vals         <- as.numeric(rna_lognorm[key, rna_neg_idx])
    gene_rna_pos[g]  <- mean(pos_vals, na.rm = TRUE)
    gene_rna_neg[g]  <- mean(neg_vals, na.rm = TRUE)
    if (length(pos_vals) > 3 && length(neg_vals) > 3)
      gene_rna_pval[g] <- wilcox.test(pos_vals, neg_vals,
                                      exact = FALSE)$p.value
  }
  rm(pk_cpm, rna_lognorm); gc()
  
  # 9. Assemble + z-score
  result <- tss_df %>%
    dplyr::select(gene_id, gene_name, display_name, chromosome, start) %>%
    dplyr::mutate(
      species    = species_label,
      group      = sex_label,
      h3k_pos    = gene_h3k_pos,
      h3k_neg    = gene_h3k_neg,
      rna_pos    = gene_rna_pos,
      rna_neg    = gene_rna_neg,
      h3k_log2fc = log2((gene_h3k_pos + 1) / (gene_h3k_neg + 1)),
      rna_log2fc = log2((gene_rna_pos + 1) / (gene_rna_neg + 1)),
      h3k_pval   = gene_h3k_pval,
      rna_pval   = gene_rna_pval
    ) %>%
    dplyr::mutate(
      z_h3k = (h3k_log2fc - mean(h3k_log2fc, na.rm = TRUE)) /
        sd(h3k_log2fc, na.rm = TRUE),
      z_rna = (rna_log2fc - mean(rna_log2fc, na.rm = TRUE)) /
        sd(rna_log2fc, na.rm = TRUE)
    )
  
  cat(sprintf("  Genes scored — H3K: %d/%d  RNA: %d/%d\n",
              sum(!is.na(result$h3k_log2fc)), n_genes,
              sum(!is.na(result$rna_log2fc)), n_genes))
  cat(sprintf("  H3K log2FC  mean=%.3f  sd=%.3f  median=%.3f\n",
              mean(result$h3k_log2fc, na.rm=TRUE),
              sd(result$h3k_log2fc,   na.rm=TRUE),
              median(result$h3k_log2fc, na.rm=TRUE)))
  cat(sprintf("  RNA log2FC  mean=%.3f  sd=%.3f  median=%.3f\n",
              mean(result$rna_log2fc, na.rm=TRUE),
              sd(result$rna_log2fc,   na.rm=TRUE),
              median(result$rna_log2fc, na.rm=TRUE)))
  
  return(result)
}

# ── Re-read barcode files (already written to disk in Block 2B) ──
bc_pl_m_pos_v <- readLines(bc_pl_m_pos)
bc_pl_m_neg_v <- readLines(bc_pl_m_neg)
bc_pl_f_pos_v <- readLines(bc_pl_f_pos)
bc_pl_f_neg_v <- readLines(bc_pl_f_neg)

cat(sprintf(
  "Pero barcodes reloaded: Male Xist+ %d / Xist- %d  |  Female Xist+ %d / Xist- %d\n",
  length(bc_pl_m_pos_v), length(bc_pl_m_neg_v),
  length(bc_pl_f_pos_v), length(bc_pl_f_neg_v)
))

# ── Rebuild tss_pl if it was also cleared ──
if (!exists("tss_pl") || nrow(tss_pl) == 0) {
  cat("Rebuilding tss_pl from GTF...\n")
  gtf_pl      <- load_gtf(gtf_path_pl, "X")
  gtf_pl$chromosome[gtf_pl$chromosome == "X"] <- "NC_051083.1"
  gtf_pl$is_x <- gtf_pl$chromosome == "NC_051083.1"
  tss_pl      <- tss_windows(gtf_pl)
  cat(sprintf("tss_pl rebuilt: %d X-linked genes\n", nrow(tss_pl)))
} else {
  cat(sprintf("tss_pl already in memory: %d genes\n", nrow(tss_pl)))
}

# ── Confirm Mus barcodes are also present for later ──
bc_mm_m_pos_v <- readLines(bc_mm_m_pos)
bc_mm_m_neg_v <- readLines(bc_mm_m_neg)
bc_mm_f_pos_v <- readLines(bc_mm_f_pos)
bc_mm_f_neg_v <- readLines(bc_mm_f_neg)

cat(sprintf(
  "Mus  barcodes reloaded: Male Xist+ %d / Xist- %d  |  Female Xist+ %d / Xist- %d\n",
  length(bc_mm_m_pos_v), length(bc_mm_m_neg_v),
  length(bc_mm_f_pos_v), length(bc_mm_f_neg_v)
))


# ── Redefine the Pero path variables (cleared from env) ──
bc_pl_m_pos <- "~/barcodes_pl_male_xist_pos.txt"
bc_pl_m_neg <- "~/barcodes_pl_male_xist_neg.txt"
bc_pl_f_pos <- "~/barcodes_pl_female_xist_pos.txt"
bc_pl_f_neg <- "~/barcodes_pl_female_xist_neg.txt"

bc_pl_m_pos_v <- readLines(bc_pl_m_pos)
bc_pl_m_neg_v <- readLines(bc_pl_m_neg)
bc_pl_f_pos_v <- readLines(bc_pl_f_pos)
bc_pl_f_neg_v <- readLines(bc_pl_f_neg)

cat(sprintf(
  "Pero barcodes reloaded: Male Xist+ %d / Xist- %d  |  Female Xist+ %d / Xist- %d\n",
  length(bc_pl_m_pos_v), length(bc_pl_m_neg_v),
  length(bc_pl_f_pos_v), length(bc_pl_f_neg_v)
))
# ============================================================
# BLOCK 3 — Run extraction for all four groups
# ============================================================

# Read barcode files (created in Block 2B)
bc_pl_m_pos_v <- readLines(bc_pl_m_pos)
bc_pl_m_neg_v <- readLines(bc_pl_m_neg)
bc_pl_f_pos_v <- readLines(bc_pl_f_pos)
bc_pl_f_neg_v <- readLines(bc_pl_f_neg)

bc_mm_m_pos_v <- readLines(bc_mm_m_pos)
bc_mm_m_neg_v <- readLines(bc_mm_m_neg)
bc_mm_f_pos_v <- readLines(bc_mm_f_pos)
bc_mm_f_neg_v <- readLines(bc_mm_f_neg)

cat(sprintf(
  "Pero barcodes: Male Xist+ %d / Xist- %d  |  Female Xist+ %d / Xist- %d\n",
  length(bc_pl_m_pos_v), length(bc_pl_m_neg_v),
  length(bc_pl_f_pos_v), length(bc_pl_f_neg_v)
))
cat(sprintf(
  "Mus  barcodes: Male Xist+ %d / Xist- %d  |  Female Xist+ %d / Xist- %d\n",
  length(bc_mm_m_pos_v), length(bc_mm_m_neg_v),
  length(bc_mm_f_pos_v), length(bc_mm_f_neg_v)
))

h3k_pero_m <- extract_paired_xist_groups(
  path_pl, tss_pl,
  bc_pl_m_pos_v, bc_pl_m_neg_v,
  "Pero", "Male"
)
saveRDS(h3k_pero_m, "~/h3k_pero_male.rds")
write.csv(h3k_pero_m, file.path(OUT_DIR, "h3k27me3_pero_male_tss2kb.csv"),
          row.names = FALSE)

h3k_pero_f <- extract_paired_xist_groups(
  path_pl, tss_pl,
  bc_pl_f_pos_v, bc_pl_f_neg_v,
  "Pero", "Female"
)
saveRDS(h3k_pero_f, "~/h3k_pero_female.rds")
write.csv(h3k_pero_f, file.path(OUT_DIR, "h3k27me3_pero_female_tss2kb.csv"),
          row.names = FALSE)

# Rebuild tss_mm with gene_id as the RNA lookup key
# (no changes needed to tss_mm itself — the function handles it)

# Re-run Mus male
h3k_mus_m <- extract_paired_xist_groups(
  path_mm, tss_mm,
  bc_mm_m_pos_v, bc_mm_m_neg_v,
  "Mus", "Male"
)
saveRDS(h3k_mus_m, "~/h3k_mus_male.rds")
write.csv(h3k_mus_m, file.path(OUT_DIR, "h3k27me3_mus_male_tss2kb.csv"),
          row.names = FALSE)

# Re-run Mus female
h3k_mus_f <- extract_paired_xist_groups(
  path_mm, tss_mm,
  bc_mm_f_pos_v, bc_mm_f_neg_v,
  "Mus", "Female"
)
saveRDS(h3k_mus_f, "~/h3k_mus_female.rds")
write.csv(h3k_mus_f, file.path(OUT_DIR, "h3k27me3_mus_female_tss2kb.csv"),
          row.names = FALSE)

# ============================================================
# BLOCK 4 — Load from disk (if re-running downstream only)
# ============================================================
if (!exists("h3k_pero_m")) h3k_pero_m <- readRDS("~/h3k_pero_male.rds")
if (!exists("h3k_pero_f")) h3k_pero_f <- readRDS("~/h3k_pero_female.rds")
if (!exists("h3k_mus_m"))  h3k_mus_m  <- readRDS("~/h3k_mus_male.rds")
if (!exists("h3k_mus_f"))  h3k_mus_f  <- readRDS("~/h3k_mus_female.rds")


# ============================================================
# BLOCK 5 — Join H3K + RNA scores to ranked DEG tables
#
# We join h3k_log2fc (raw, within-species) for within-species
# plots and z_h3k / z_rna for cross-species comparisons.
# ============================================================
cat("\nJoining H3K27me3 + RNA scores to DEG tables...\n")

join_paired <- function(ranked_df, scored_df) {
  ranked_df %>%
    dplyr::left_join(
      scored_df %>%
        dplyr::select(gene_name,
                      h3k_pos, h3k_neg,
                      rna_pos, rna_neg,
                      h3k_log2fc, rna_log2fc,
                      z_h3k,     z_rna,
                      h3k_pval,  rna_pval),
      by = "gene_name"
    )
}

ranked_mus_m  <- join_paired(ranked_mus_m,  h3k_mus_m)
ranked_mus_f  <- join_paired(ranked_mus_f,  h3k_mus_f)
ranked_pero_m <- join_paired(ranked_pero_m, h3k_pero_m)
ranked_pero_f <- join_paired(ranked_pero_f, h3k_pero_f)

cat("H3K27me3 join coverage:\n")
for (nm in c("ranked_mus_m","ranked_mus_f","ranked_pero_m","ranked_pero_f")) {
  df <- get(nm)
  cat(sprintf("  %-20s %d / %d genes with H3K log2FC\n",
              nm, sum(!is.na(df$h3k_log2fc)), nrow(df)))
}


# ============================================================
# BLOCK 6 — Ranked lollipop plots
# ============================================================
cat("\nPlotting ranked lollipop plots...\n")

plot_ranked_lollipop <- function(df, title, top_n = 50) {
  
  df <- df %>%
    
    # --------------------------------------------------------
  # Ensure clean numeric values
  # --------------------------------------------------------
  mutate(
    h3k = h3k_log2fc,
    rna = rna_log2fc
  ) %>%
    
    # --------------------------------------------------------
  # Rank by MOST SILENCED in Xist+ (most negative RNA log2FC)
  # --------------------------------------------------------
  arrange(rna) %>%
    mutate(rank = row_number()) %>%
    slice_head(n = top_n)
  
  ggplot(df, aes(x = rank, y = rna)) +
    
    # lollipop stems
    geom_segment(aes(xend = rank, y = 0, yend = rna),
                 color = "grey80", linewidth = 0.4) +
    
    # points
    geom_point(
      aes(
        size = abs(h3k) + 0.1,
        fill = h3k
      ),
      shape = 21,
      color = "grey30",
      alpha = 0.9
    ) +
    
    # gene labels
    geom_text_repel(
      aes(label = gene_name),
      size               = 2.6,
      fontface           = "italic",
      color              = "grey20",
      max.overlaps       = Inf,
      box.padding        = 0.3,
      point.padding      = 0.2,
      segment.size       = 0.25,
      segment.color      = "grey60",
      segment.alpha      = 0.6,
      min.segment.length = 0.1,
      force              = 2
    ) +
    
    # H3K color scale
    scale_fill_gradient2(
      low      = "steelblue",
      mid      = "white",
      high     = "#b2182b",
      midpoint = 0,
      name     = "H3K27me3\nlog2FC",
      na.value = "grey70"
    ) +
    
    scale_size_continuous(range = c(1.5, 7), guide = "none") +
    
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    
    labs(
      title = title,
      x     = "Rank (1 = most silenced in Xist+)",
      y     = expression("RNA log"[2]*"FC (Xist+ / Xist-)") 
    ) +
    
    pub_theme +
    
    theme(
      legend.position = "right",
      plot.title = element_text(size = BASE, face = "bold",
                                hjust = 0.5, color = "black")
    )
}

ranked_mus_m <- paired_all_fixed %>%
  filter(species == "Mus", group == "Male")

ranked_mus_f <- paired_all_fixed %>%
  filter(species == "Mus", group == "Female")

ranked_pero_m <- paired_all_fixed %>%
  filter(species == "Pero", group == "Male")

ranked_pero_f <- paired_all_fixed %>%
  filter(species == "Pero", group == "Female")
fig_lollipop <- (p_lol_mm_m | p_lol_mm_f) / (p_lol_pl_m | p_lol_pl_f) +
  plot_annotation(
    title      = "Silenced X-linked Genes Ranked by log2FC",
    subtitle   = "Dot size + color = H3K27me3 log2FC at TSS ± 2kb (CPM, within-species)",
    tag_levels = "a",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
      plot.tag      = element_text(size = BASE + 2, face = "bold")
    )
  )

ggsave(file.path(OUT_DIR, "ranked_lollipop_h3k27me3.pdf"),
       fig_lollipop, width = 20, height = 16, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "ranked_lollipop_h3k27me3.png"),
       fig_lollipop, width = 20, height = 16, units = "in", dpi = 300)
cat("Saved ranked_lollipop_h3k27me3.pdf/.png\n")


# ============================================================
# BLOCK 7 — Scatter: RNA silencing vs H3K27me3
#   Uses z-scored values for cross-species panels
# ============================================================
cat("\nPlotting FC vs H3K27me3 scatter...\n")

all_scatter <- bind_rows(
  ranked_mus_m  %>% mutate(group = "Mus Male"),
  ranked_mus_f  %>% mutate(group = "Mus Female"),
  ranked_pero_m %>% mutate(group = "Pero Male"),
  ranked_pero_f %>% mutate(group = "Pero Female")
) %>%
  filter(!is.na(h3k_log2fc)) %>%
  mutate(group = factor(group,
                        levels = c("Mus Male","Mus Female",
                                   "Pero Male","Pero Female")))

all_scatter <- all_scatter %>%
  mutate(
    category = case_when(
      avg_log2FC < 0  & h3k_log2fc > 0  ~ "Silenced + H3K27me3 gain",
      avg_log2FC < 0  & h3k_log2fc <= 0 ~ "Silenced only",
      avg_log2FC >= 0 & h3k_log2fc > 0  ~ "H3K27me3 gain only",
      TRUE                               ~ "Other"
    )
  )

r_df <- all_scatter %>%
  group_by(group) %>%
  summarise(r = cor(avg_log2FC, h3k_log2fc, use = "complete.obs"),
            .groups = "drop") %>%
  mutate(label = sprintf("r = %.2f", r),
         x = Inf, y = Inf)

fig_scatter <- ggplot(all_scatter,
                      aes(x = avg_log2FC, y = h3k_log2fc, color = category)) +
  geom_point(data = ~ filter(.x, category == "Other"),
             size = 1.2, alpha = 0.4) +
  geom_point(data = ~ filter(.x, category != "Other"),
             size = 1.8, alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_text_repel(
    aes(label = gene_name),
    size               = 2.4,
    fontface           = "italic",
    segment.size       = 0.25,
    segment.alpha      = 0.5,
    max.overlaps       = Inf,
    box.padding        = 0.3,
    point.padding      = 0.15,
    min.segment.length = 0.1,
    force              = 2
  ) +
  geom_text(data        = r_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE,
            hjust = 1.15, vjust = 1.5,
            size = BASE * 0.28, color = "grey30") +
  scale_color_manual(values = cat_colors, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  scale_y_continuous(expand = expansion(mult = 0.08),
                     labels = number_format(accuracy = 0.1)) +
  facet_wrap(~ group, nrow = 2, scales = "free") +
  labs(
    x = expression(italic("Xist+")*" vs "*italic("Xist-")*" RNA log"[2]*"FC"),
    y = expression(italic("Xist+")*" vs "*italic("Xist-")*" H3K27me3 log"[2]*"FC (CPM)")
  ) +
  pub_theme +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = BASE - 2),
    strip.text      = element_text(face = "bold", size = BASE),
    panel.border    = element_rect(color = "grey80", fill = NA, linewidth = 0.4)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

ggsave(file.path(OUT_DIR, "scatter_silencing_h3k27me3.pdf"),
       fig_scatter, width = 12, height = 11, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "scatter_silencing_h3k27me3.png"),
       fig_scatter, width = 12, height = 11, units = "in", dpi = 300)
cat("Saved scatter_silencing_h3k27me3.pdf/.png\n")

labelled_genes <- all_scatter %>%
  filter(category == "Silenced + H3K27me3 gain") %>%
  dplyr::select(group, gene_name, avg_log2FC, h3k_log2fc, rank) %>%
  arrange(group, avg_log2FC)
write.csv(labelled_genes,
          file.path(OUT_DIR, "silenced_h3k27me3_high_genes.csv"),
          row.names = FALSE)
cat(sprintf("Silenced + H3K27me3 gain genes: %d total\n", nrow(labelled_genes)))
print(labelled_genes %>% count(group))


# ============================================================
# BLOCK 8 — Cross-species ortholog comparison
# ============================================================
cat("\nBuilding cross-species ortholog comparison...\n")

mm_m_join <- ranked_mus_m  %>%
  dplyr::select(gene_name,
                fc_mus_m    = avg_log2FC,
                h3k_mus_m   = h3k_log2fc,
                z_h3k_mus_m = z_h3k,
                rank_mus_m  = rank)
mm_f_join <- ranked_mus_f  %>%
  dplyr::select(gene_name, fc_mus_f = avg_log2FC, rank_mus_f = rank)
pl_m_join <- ranked_pero_m %>%
  dplyr::select(gene_name,
                fc_pero_m    = avg_log2FC,
                h3k_pero_m   = h3k_log2fc,
                z_h3k_pero_m = z_h3k,
                rank_pero_m  = rank)
pl_f_join <- ranked_pero_f %>%
  dplyr::select(gene_name, fc_pero_f = avg_log2FC, rank_pero_f = rank)

orth_all <- orthologs %>%
  left_join(mm_m_join, by = c("shared_name" = "gene_name")) %>%
  left_join(mm_f_join, by = c("shared_name" = "gene_name")) %>%
  left_join(pl_m_join, by = c("shared_name" = "gene_name")) %>%
  left_join(pl_f_join, by = c("shared_name" = "gene_name")) %>%
  mutate(
    sil_mus_m  = !is.na(fc_mus_m)  & fc_mus_m  < 0,
    sil_mus_f  = !is.na(fc_mus_f)  & fc_mus_f  < 0,
    sil_pero_m = !is.na(fc_pero_m) & fc_pero_m < 0,
    sil_pero_f = !is.na(fc_pero_f) & fc_pero_f < 0,
    n_groups_silenced      = sil_mus_m + sil_mus_f + sil_pero_m + sil_pero_f,
    conserved_both_species = (sil_mus_m & sil_pero_m) | (sil_mus_f & sil_pero_f)
  ) %>%
  arrange(desc(n_groups_silenced), fc_mus_m)

write.csv(orth_all,
          file.path(OUT_DIR, "ortholog_silencing_comparison.csv"),
          row.names = FALSE)
cat(sprintf("Conserved silenced (both species): %d genes\n",
            sum(orth_all$conserved_both_species, na.rm = TRUE)))

plot_orth_scatter <- function(df, x_col, y_col, x_lab, y_lab) {
  df <- df %>% filter(!is.na(.data[[x_col]]), !is.na(.data[[y_col]]))
  r  <- cor(df[[x_col]], df[[y_col]], use = "complete.obs")
  df <- df %>%
    mutate(both_sil = .data[[x_col]] < 0 & .data[[y_col]] < 0)
  
  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_point(aes(color = both_sil), size = 2, alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, color = "#B5294E",
                fill = "#E8A0A0", linewidth = 0.8, alpha = 0.2) +
    geom_text_repel(
      aes(label = shared_name),
      size               = 2.5,
      fontface           = "italic",
      color              = "grey20",
      max.overlaps       = Inf,
      box.padding        = 0.3,
      point.padding      = 0.15,
      segment.size       = 0.25,
      segment.color      = "grey60",
      segment.alpha      = 0.5,
      min.segment.length = 0.1,
      force              = 2
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#B5294E"),
                       guide  = "none") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = sprintf("r = %.2f", r), size = 4, color = "#B5294E") +
    labs(x = x_lab, y = y_lab) +
    pub_theme
}

p_orth_male    <- plot_orth_scatter(orth_all, "fc_mus_m",  "fc_pero_m",
                                    "Mus Male log2FC",   "Pero Male log2FC")
p_orth_female  <- plot_orth_scatter(orth_all, "fc_mus_f",  "fc_pero_f",
                                    "Mus Female log2FC", "Pero Female log2FC")
p_orth_mf_mus  <- plot_orth_scatter(orth_all, "fc_mus_m",  "fc_mus_f",
                                    "Mus Male log2FC",   "Mus Female log2FC")
p_orth_mf_pero <- plot_orth_scatter(orth_all, "fc_pero_m", "fc_pero_f",
                                    "Pero Male log2FC",  "Pero Female log2FC")

fig_ortholog <- (p_orth_male | p_orth_female) / (p_orth_mf_mus | p_orth_mf_pero) +
  plot_annotation(
    title      = "Cross-Species X-linked Gene Silencing (Ortholog-Matched)",
    subtitle   = "Red = silenced in both species",
    tag_levels = "a",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
      plot.tag      = element_text(size = BASE + 2, face = "bold")
    )
  )

ggsave(file.path(OUT_DIR, "ortholog_crossspecies_scatter.pdf"),
       fig_ortholog, width = 16, height = 14, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "ortholog_crossspecies_scatter.png"),
       fig_ortholog, width = 16, height = 14, units = "in", dpi = 300)
cat("Saved ortholog_crossspecies_scatter.pdf/.png\n")


# ============================================================
# BLOCK 9 — Venn diagrams
# ============================================================
cat("\nBuilding Venn diagrams...\n")

sil_sets <- list(
  "Mus Male"    = ranked_mus_m$gene_name[ ranked_mus_m$avg_log2FC   <= -0.58],
  "Mus Female"  = ranked_mus_f$gene_name[ ranked_mus_f$avg_log2FC   <= -0.58],
  "Pero Male"   = ranked_pero_m$gene_name[ranked_pero_m$avg_log2FC  <= -0.58],
  "Pero Female" = ranked_pero_f$gene_name[ranked_pero_f$avg_log2FC  <= -0.58]
)

make_venn2 <- function(set_list, title, hi_color) {
  ggVennDiagram(set_list, label_alpha = 0,
                edge_size = 0.7, label_size = 3.5) +
    scale_fill_gradient(low = "#F5F5F5", high = hi_color) +
    scale_color_manual(values = rep("grey40", length(set_list))) +
    labs(title = title) +
    theme_void(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5, size = 11,
                                     margin = margin(b = 6)),
      plot.margin     = margin(30, 30, 30, 30),
      legend.position = "none"
    )
}

p_v1 <- make_venn2(list("Mus Male"    = sil_sets[["Mus Male"]],
                        "Pero Male"   = sil_sets[["Pero Male"]]),
                   "Silenced in Males\nMus vs Pero",   "#2171b5")
p_v2 <- make_venn2(list("Mus Female"  = sil_sets[["Mus Female"]],
                        "Pero Female" = sil_sets[["Pero Female"]]),
                   "Silenced in Females\nMus vs Pero", "#cb181d")
p_v3 <- make_venn2(list("Mus Male"    = sil_sets[["Mus Male"]],
                        "Mus Female"  = sil_sets[["Mus Female"]]),
                   "Male vs Female\nWithin Mus",       "#238b45")
p_v4 <- make_venn2(list("Pero Male"   = sil_sets[["Pero Male"]],
                        "Pero Female" = sil_sets[["Pero Female"]]),
                   "Male vs Female\nWithin Pero",      "#6a51a3")

p_v_4way <- ggVennDiagram(sil_sets, label_alpha = 0,
                          edge_size = 0.7, label_size = 3) +
  scale_fill_gradient(low = "#F5F5F5", high = "#B5294E") +
  labs(title = "All Groups") +
  theme_void(base_size = 11) +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5, size = 12),
        legend.position = "none")

fig_venns <- (((p_v1 | p_v2) / (p_v3 | p_v4)) | wrap_elements(p_v_4way)) +
  plot_layout(widths = c(2, 1)) +
  plot_annotation(
    title      = "Silenced X-linked DEGs: Overlap Across Groups",
    subtitle   = "avg_log2FC \u2264 -0.58 (>1.5\u00d7 silenced)",
    tag_levels = "a",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
      plot.tag      = element_text(size = BASE + 2, face = "bold")
    )
  )

ggsave(file.path(OUT_DIR, "venn_silenced_xlinked.pdf"),
       fig_venns, width = 18, height = 12, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "venn_silenced_xlinked.png"),
       fig_venns, width = 18, height = 12, units = "in", dpi = 300)
cat("Saved venn_silenced_xlinked.pdf/.png\n")

cat("\n── Overlap summaries ──\n")
for (a in names(sil_sets)) {
  for (b in names(sil_sets)) {
    if (a >= b) next
    ov <- intersect(sil_sets[[a]], sil_sets[[b]])
    cat(sprintf("\n%s | %s\n  Overlap (%d): %s\n",
                a, b, length(ov),
                if (length(ov) > 0) paste(sort(ov), collapse = ", ") else "none"))
  }
}


# ============================================================
# BLOCK 10 — Focal gene spotlight
# ============================================================
cat("\nPlotting focal gene spotlight...\n")

focal_genes <- c("Awat2", "Bex3", "Armcx4", "Lancl3", "Drp2")

pull_focal <- function(df, sex_label) {
  df %>%
    filter(gene_name %in% focal_genes) %>%
    mutate(
      sex       = sex_label,
      sig       = sapply(p_val_adj, fmt_p),
      gene_name = factor(gene_name, levels = focal_genes)
    ) %>%
    dplyr::select(gene_name, sex, avg_log2FC, pct.1, pct.2, p_val_adj, sig)
}

focal_pero_m <- pull_focal(deg_pero_m, "Male")
focal_pero_f <- pull_focal(deg_pero_f, "Female")

missing_m <- setdiff(focal_genes, as.character(focal_pero_m$gene_name))
missing_f <- setdiff(focal_genes, as.character(focal_pero_f$gene_name))
if (length(missing_m) > 0)
  warning("Missing from pero_male DEG: ", paste(missing_m, collapse = ", "))
if (length(missing_f) > 0)
  warning("Missing from pero_female DEG: ", paste(missing_f, collapse = ", "))

focal <- bind_rows(focal_pero_m, focal_pero_f) %>%
  mutate(sex = factor(sex, levels = c("Male", "Female")))

cat("\nFocal gene stats:\n")
focal %>%
  dplyr::select(gene_name, sex, avg_log2FC, pct.1, pct.2, p_val_adj, sig) %>%
  arrange(sex, avg_log2FC) %>%
  print(n = Inf)

dot_data <- bind_rows(
  focal %>% transmute(gene_name, sex,
                      group   = paste(as.character(sex), "Xist+"),
                      pct_exp = pct.1 * 100,
                      log2fc  = avg_log2FC,
                      sig     = sig),
  focal %>% transmute(gene_name, sex,
                      group   = paste(as.character(sex), "Xist-"),
                      pct_exp = pct.2 * 100,
                      log2fc  = 0,
                      sig     = "")
) %>%
  mutate(
    group    = factor(group,
                      levels = c("Male Xist+","Male Xist-",
                                 "Female Xist+","Female Xist-")),
    xist_pos = grepl("Xist\\+", group)
  )

fc_lim <- ceiling(max(abs(focal$avg_log2FC), na.rm = TRUE) * 10) / 10

pA <- ggplot(dot_data,
             aes(x = group, y = gene_name,
                 size = pct_exp, fill = log2fc)) +
  geom_point(shape = 21, color = "grey25", stroke = 0.35, alpha = 0.92) +
  geom_text(
    data    = dot_data %>% filter(sig != "" & xist_pos),
    aes(label = sig),
    size    = 3.8, vjust = -1.2, fontface = "bold",
    color   = "grey10", show.legend = FALSE
  ) +
  geom_vline(xintercept = 2.5,
             linetype = "dashed", color = "grey55", linewidth = 0.4) +
  annotate("text", x = 1.5, y = length(focal_genes) + 0.85,
           label = "Male",   fontface = "bold", size = 4.2, color = "#2171b5") +
  annotate("text", x = 3.5, y = length(focal_genes) + 0.85,
           label = "Female", fontface = "bold", size = 4.2, color = "#cb181d") +
  scale_size_continuous(
    range  = c(2, 11),
    breaks = c(10, 25, 50, 75),
    name   = "% Cells\nExpressing"
  ) +
  scale_fill_gradient2(
    low      = "#2166AC", mid = "grey97", high = "#B5294E",
    midpoint = 0, limits = c(-fc_lim, fc_lim),
    name     = expression("log"[2]*"FC\n(Xist"^"+"~"/ Xist"^"-"*")")
  ) +
  scale_x_discrete(labels = c("Xist+", "Xist-", "Xist+", "Xist-")) +
  scale_y_discrete(limits = rev(focal_genes)) +
  labs(x = NULL, y = NULL,
       title = "Pero: Expression by Sex and Xist Status") +
  pub_theme +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = BASE - 2),
    axis.text.y      = element_text(face = "italic"),
    panel.grid.major = element_line(color = "grey91", linewidth = 0.3)
  ) +
  guides(
    fill = guide_colorbar(barwidth = 0.7, barheight = 5,
                          ticks = FALSE, title.hjust = 0.5),
    size = guide_legend(override.aes = list(fill = "grey60", color = "grey25"))
  )

gene_order <- focal_pero_m %>%
  arrange(avg_log2FC) %>%
  pull(gene_name) %>%
  as.character()

focal_plot <- focal %>%
  mutate(
    gene_name = factor(gene_name, levels = gene_order),
    hjust_val = ifelse(avg_log2FC < 0, 1, 0),
    nudge_val = ifelse(avg_log2FC < 0, -0.07, 0.07)
  )

pB <- ggplot(focal_plot,
             aes(x = avg_log2FC, y = gene_name, fill = sex)) +
  geom_col(position = position_dodge(width = 0.65),
           width = 0.55, color = NA, alpha = 0.88) +
  geom_text(
    aes(x     = avg_log2FC + nudge_val,
        label = sig,
        group = sex),
    position  = position_dodge(width = 0.65),
    hjust     = focal_plot$hjust_val,
    size      = 3.8, fontface = "bold",
    color     = "grey15", show.legend = FALSE
  ) +
  geom_vline(xintercept = 0,     color = "black",   linewidth = 0.45) +
  geom_vline(xintercept = -0.58, linetype = "dashed",
             color = "grey50",  linewidth = 0.35) +
  annotate("text", x = -0.61, y = 0.45,
           label = "1.5\u00d7 threshold", angle = 90,
           hjust = 0, size = 2.8, color = "grey45") +
  scale_fill_manual(
    values = c("Male" = "#2171b5", "Female" = "#cb181d"),
    name   = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = 0.16)) +
  labs(
    x     = expression("avg log"[2]*"FC  (Xist"^"+"~"vs Xist"^"-"*")"),
    y     = NULL,
    title = "Silencing Magnitude by Sex"
  ) +
  pub_theme +
  theme(
    legend.position    = "top",
    axis.text.y        = element_text(face = "italic"),   # fixed broken syntax
    panel.grid.major.x = element_line(color = "grey91", linewidth = 0.3)
  )

fig_focal <- (pA | pB) +
  plot_annotation(
    title      = "X-linked Genes Silenced in Pero Xist\u207a Males",
    subtitle   = "Col4a5 \u00b7 Lonrf3 \u00b7 Mct1 \u00b7 Reps2 \u00b7 Acsl4 \u00b7 Dcx  \u2014  bulk + snRNA validated",
    tag_levels = "a",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
      plot.tag      = element_text(size = BASE + 2, face = "bold")
    )
  ) +
  plot_layout(widths = c(1.6, 1))

ggsave(file.path(OUT_DIR, "focal_genes_pero_male_silenced.pdf"),
       fig_focal, width = 13, height = 6, units = "in", device = cairo_pdf)
ggsave(file.path(OUT_DIR, "focal_genes_pero_male_silenced.png"),
       fig_focal, width = 13, height = 6, units = "in", dpi = 300)
cat("Saved focal_genes_pero_male_silenced.pdf/.png\n")


# ============================================================
# BLOCK 11 — Selected gene absolute barplot (within-species)
#
# Uses raw CPM means (h3k_pos/h3k_neg) and log-norm means
# (rna_pos/rna_neg) — no z-scoring — so y-axis is interpretable
# in native units for each species separately.
#
# For cross-species comparisons use z_h3k / z_rna columns.
# ============================================================

# Build a combined long table with both species and both sexes
paired_all <- bind_rows(
  h3k_pero_m %>% mutate(dataset = "Pero"),
  h3k_pero_f %>% mutate(dataset = "Pero"),
  h3k_mus_m  %>% mutate(dataset = "Mus"),
  h3k_mus_f  %>% mutate(dataset = "Mus")
)
write.csv(paired_all, file.path(OUT_DIR, "paired_all.csv"), row.names = FALSE)
#' Plot absolute Xist+/Xist- RNA and H3K27me3 for selected genes
#'
#' @param df       Combined paired_all table
#' @param genes    Character vector of gene display_names to show
#' @param datasets Which datasets to include: "Pero", "Mus", or both
#' @param groups   Which sex groups: "Male", "Female", or both
#'
#' Within-species: set datasets to one species, scales are free_x
#' Cross-species:  set datasets to both; values are z-scored (z_h3k/z_rna)
plot_selected_genes_absolute <- function(df,
                                         genes,
                                         datasets   = c("Pero", "Mus"),
                                         groups     = c("Female", "Male"),
                                         use_zscore = FALSE,
                                         xlim_rna   = NULL,   # <-- ADD
                                         xlim_h3k   = NULL) { # <-- ADD {
  library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
  
  value_h3k_pos <- if (use_zscore) "z_h3k" else "h3k_pos"
  value_h3k_neg <- if (use_zscore) NA       else "h3k_neg"
  value_rna_pos <- if (use_zscore) "z_rna"  else "rna_pos"
  value_rna_neg <- if (use_zscore) NA        else "rna_neg"
  x_lab_h3k     <- if (use_zscore) "H3K27me3 z-score (Xist+ log2FC)" else "H3K27me3 mean CPM"
  x_lab_rna     <- if (use_zscore) "RNA z-score (Xist+ log2FC)"      else "RNA log-norm mean"
  
  df_plot <- df %>%
    dplyr::filter(display_name %in% genes,
                  dataset %in% datasets,
                  group   %in% groups) %>%
    # With this — keeps rows that have at least RNA data
    dplyr::filter(!is.na(rna_pos)) %>%   
    dplyr::mutate(display_name = factor(display_name, levels = rev(genes)),
                  facet_label  = paste(dataset, group))
  
  if (nrow(df_plot) == 0) {
    warning("No rows after filtering — check gene names, datasets, groups")
    return(invisible(NULL))
  }
  
  if (!use_zscore) {
    # absolute mode: show Xist+ and Xist- as separate bars
    df_rna <- df_plot %>%
      dplyr::select(display_name, dataset, group, facet_label, rna_pos, rna_neg) %>%
      pivot_longer(c(rna_pos, rna_neg),
                   names_to  = "type", values_to = "value") %>%
      mutate(status = ifelse(type == "rna_pos", "Xist+", "Xist-"))
    
    df_h3k <- df_plot %>%
      dplyr::select(display_name, dataset, group, facet_label, h3k_pos, h3k_neg) %>%
      pivot_longer(c(h3k_pos, h3k_neg),
                   names_to  = "type", values_to = "value") %>%
      mutate(status = ifelse(type == "h3k_pos", "Xist+", "Xist-"))
    
    p_rna <- ggplot(df_rna, aes(x = value, y = display_name, fill = status)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      facet_grid(. ~ facet_label, scales = "free_x") +
      scale_fill_manual(values = c("Xist+" = "#4CAF50", "Xist-" = "#E64B35"), name = NULL) +
      labs(x = x_lab_rna, y = NULL, title = "RNA expression") +
      { if (!is.null(xlim_rna)) coord_cartesian(xlim = xlim_rna) } +  # <-- ADD
      pub_theme +
      theme(axis.text.y = element_text(face = "italic"),
            legend.position = "top",
            strip.text = element_text(face = "bold"))
    
    p_h3k <- ggplot(df_h3k, aes(x = value, y = display_name, fill = status)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      facet_grid(. ~ facet_label, scales = "free_x") +
      scale_fill_manual(values = c("Xist+" = "#4CAF50", "Xist-" = "#E64B35"), name = NULL) +
      labs(x = x_lab_h3k, y = NULL, title = "H3K27me3 signal") +
      { if (!is.null(xlim_h3k)) coord_cartesian(xlim = xlim_h3k) } +  # <-- ADD
      pub_theme +
      theme(axis.text.y = element_text(face = "italic"),
            legend.position = "none",
            strip.text = element_text(face = "bold"))
    
  } else {
    # z-score mode: single bar = log2FC z-score, diverging colour
    df_rna <- df_plot %>%
      dplyr::select(display_name, dataset, group, facet_label, z_rna)
    
    df_h3k <- df_plot %>%
      dplyr::select(display_name, dataset, group, facet_label, z_h3k)
    
    p_rna <- ggplot(df_rna,
                    aes(x = z_rna, y = display_name, fill = z_rna)) +
      geom_col(width = 0.6) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
      facet_grid(. ~ facet_label, scales = "free_x") +
      scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B5294E",
                           midpoint = 0, guide = "none") +
      labs(x = x_lab_rna, y = NULL, title = "RNA (z-scored log2FC)") +
      pub_theme +
      theme(axis.text.y  = element_text(face = "italic"),
            strip.text   = element_text(face = "bold"))
    
    p_h3k <- ggplot(df_h3k,
                    aes(x = z_h3k, y = display_name, fill = z_h3k)) +
      geom_col(width = 0.6) +
      geom_vline(xintercept = 0, color = "black", linewidth = 0.4) +
      facet_grid(. ~ facet_label, scales = "free_x") +
      scale_fill_gradient2(low = "#2166AC", mid = "grey90", high = "#B5294E",
                           midpoint = 0, guide = "none") +
      labs(x = x_lab_h3k, y = NULL, title = "H3K27me3 (z-scored log2FC)") +
      pub_theme +
      theme(axis.text.y  = element_text(face = "italic"),
            strip.text   = element_text(face = "bold"))
  }
  
  p_rna / p_h3k
}
# Compute shared x-axis limits across both species
df_lim <- paired_all %>%
  filter(display_name %in% genes_of_interest,
         dataset %in% c("Pero", "Mus"),
         group   %in% c("Female", "Male"),
         !is.na(h3k_pos), !is.na(h3k_neg))

shared_xlim_rna <- c(0, max(c(df_lim$rna_pos, df_lim$rna_neg), na.rm = TRUE) * 1.05)
shared_xlim_h3k <- c(0, max(c(df_lim$h3k_pos, df_lim$h3k_neg), na.rm = TRUE) * 1.05)

# Pero
fig_pero_abs <- plot_selected_genes_absolute(
  paired_all, genes_of_interest,
  datasets   = "Pero",
  groups     = c("Female", "Male"),
  use_zscore = FALSE,
  xlim_rna   = shared_xlim_rna,
  xlim_h3k   = shared_xlim_h3k
)

# Mus
fig_mus_abs <- plot_selected_genes_absolute(
  paired_all, genes_of_interest,
  datasets   = "Mus",
  groups     = c("Female", "Male"),
  use_zscore = FALSE,
  xlim_rna   = shared_xlim_rna,
  xlim_h3k   = shared_xlim_h3k
)

fig_pero_abs | fig_mus_abs

# ── Example usage ────────────────────────────────────────────

genes_of_interest <- c("Awat2", "Bex3", "Armcx4", "Lancl3", "Drp2","")

# Within-species absolute (Pero only)
fig_pero_abs <- plot_selected_genes_absolute(
  paired_all, genes_of_interest,
  datasets   = "Pero",
  groups     = c("Female", "Male"),
  use_zscore = FALSE
)
ggsave(file.path(OUT_DIR, "focal_pero_absolute_barplot.pdf"),
       fig_pero_abs, width = 10, height = 7, units = "in", device = cairo_pdf)

# Within-species absolute (Mus only)
fig_mus_abs <- plot_selected_genes_absolute(
  paired_all, genes_of_interest,
  datasets   = "Mus",
  groups     = c("Female", "Male"),
  use_zscore = FALSE
)
ggsave(file.path(OUT_DIR, "focal_mus_absolute_barplot.pdf"),
       fig_mus_abs, width = 10, height = 7, units = "in", device = cairo_pdf)
fig_pero_abs | fig_mus_abs
# Cross-species z-scored (both species together)
fig_cross_z <- plot_selected_genes_absolute(
  paired_all, genes_of_interest,
  datasets   = c("Pero", "Mus"),
  groups     = c("Female", "Male"),
  use_zscore = TRUE
)
ggsave(file.path(OUT_DIR, "focal_cross_species_zscore_barplot.pdf"),
       fig_cross_z, width = 14, height = 7, units = "in", device = cairo_pdf)

cat("Saved focal gene barplots\n")


# ============================================================
# BLOCK 12 — Bulk + snPairedTag integration
# ============================================================
cat("\nRunning bulk vs snPairedTag integration...\n")

# Set A: all bulk genes (already all pero male silenced)
# (assumes `bulk` data frame exists with gene_name column)
set_bulk <- bulk$gene_name
cat(sprintf("Set A — Bulk Pero male silenced: %d genes\n", length(set_bulk)))

set_sn_pero_m_rna <- deg_pero_m$gene_name[deg_pero_m$avg_log2FC < 0]
cat(sprintf("Set B — snPairedTag Pero male RNA silenced: %d genes\n",
            length(set_sn_pero_m_rna)))

set_sn_pero_m_h3k <- h3k_pero_m$gene_name[
  !is.na(h3k_pero_m$h3k_log2fc) & h3k_pero_m$h3k_log2fc > 0
]
cat(sprintf("Set C — snPairedTag Pero male H3K27me3 gain: %d genes\n",
            length(set_sn_pero_m_h3k)))

set_sn_pero_f_rna <- deg_pero_f$gene_name[deg_pero_f$avg_log2FC < 0]
set_sn_mus_m_rna  <- deg_mus_m$gene_name[deg_mus_m$avg_log2FC < 0]

core_genes  <- Reduce(intersect,
                      list(set_bulk, set_sn_pero_m_rna, set_sn_pero_m_h3k))
bulk_sn_rna <- intersect(set_bulk, set_sn_pero_m_rna)
bulk_sn_h3k <- intersect(set_bulk, set_sn_pero_m_h3k)

cat(sprintf(
  "\nOverlap summary:\n  Bulk ∩ sn RNA (FC<0):   %d\n  Bulk ∩ sn H3K (FC>0):   %d\n  CORE (all three):        %d\n",
  length(bulk_sn_rna), length(bulk_sn_h3k), length(core_genes)
))
cat("Core genes:", paste(sort(core_genes), collapse = ", "), "\n")

integrated <- bulk %>%
  left_join(deg_pero_m %>% dplyr::select(gene_name, sn_pero_m_rna_fc = avg_log2FC), by = "gene_name") %>%
  left_join(deg_pero_f %>% dplyr::select(gene_name, sn_pero_f_rna_fc = avg_log2FC), by = "gene_name") %>%
  left_join(deg_mus_m  %>% dplyr::select(gene_name, sn_mus_m_rna_fc  = avg_log2FC), by = "gene_name") %>%
  left_join(deg_mus_f  %>% dplyr::select(gene_name, sn_mus_f_rna_fc  = avg_log2FC), by = "gene_name") %>%
  left_join(h3k_pero_m %>% dplyr::select(gene_name, h3k_sn_pero_m = h3k_log2fc),   by = "gene_name") %>%
  left_join(h3k_pero_f %>% dplyr::select(gene_name, h3k_sn_pero_f = h3k_log2fc),   by = "gene_name") %>%
  left_join(h3k_mus_m  %>% dplyr::select(gene_name, h3k_sn_mus_m  = h3k_log2fc),   by = "gene_name") %>%
  left_join(h3k_mus_f  %>% dplyr::select(gene_name, h3k_sn_mus_f  = h3k_log2fc),   by = "gene_name") %>%
  mutate(
    in_sn_pero_m_rna = !is.na(sn_pero_m_rna_fc) & sn_pero_m_rna_fc < 0,
    in_sn_pero_m_h3k = !is.na(h3k_sn_pero_m)   & h3k_sn_pero_m   > 0,
    in_sn_pero_f_rna = !is.na(sn_pero_f_rna_fc) & sn_pero_f_rna_fc < 0,
    in_sn_mus_m_rna  = !is.na(sn_mus_m_rna_fc)  & sn_mus_m_rna_fc  < 0,
    is_core          = gene_name %in% core_genes,
    n_sn_evidence    = as.integer(in_sn_pero_m_rna) + as.integer(in_sn_pero_m_h3k)
  ) %>%
  arrange(desc(is_core), desc(n_sn_evidence))

write.csv(integrated,
          file.path(OUT_DIR, "integrated_bulk_sn_h3k_pero_male.csv"),
          row.names = FALSE)
cat(sprintf("Integrated table saved: %d genes\n", nrow(integrated)))

cat("\n=== ALL DONE ===\n")
cat("Outputs in", OUT_DIR, "\n")