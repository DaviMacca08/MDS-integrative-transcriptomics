# =========================================================
#      Cross-dataset comparison of SVA-corrected DEGs
# =========================================================
#
# Description:
# This script compares the number of significantly differentially
# expressed genes (DEGs) identified after surrogate variable
# analysis (SVA) correction across three independent MDS
# transcriptomic datasets:
#
#   - GSE19429 (Microarray)
#   - GSE58831 (Microarray)
#   - GSE114922 (RNA-seq)
#
# Differential expression results obtained using the Buja &
# Eyuboglu (BE) and Leek methods are summarized and visualized
# as grouped bar plots, reporting the number of upregulated and
# downregulated genes in MDS.
#
# Significant DEGs are defined as:
#   - Adjusted P-value < 0.05
#   - |log2 Fold Change| > 1
#
# Notes:
# - The Leek method is not available for GSE114922 due to SVA
#   convergence failure.
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("SetupEnvironment/00_paths.R")
source("SetupEnvironment/01_environment.R")
source("SetupEnvironment/02_seed.R")
source("SetupEnvironment/03_helper_functions.R")

message("=== Cross-dataset comparison: GSE19429 vs GSE58831 vs GSE114922 ===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Reading SVA-corrected differential expression results...")

deg_gse19429_sva_be <- read.csv(file.path(paths$tables_sva,  "DEGs_SVA_BE_GSE19429.csv"), sep = ",")
deg_gse19429_sva_leek <- read.csv(file.path(paths$tables_sva,  "DEGs_SVA_LEEK_GSE19429.csv"), sep = ",")

deg_gse58831_sva_be <- read.csv(file.path(paths$tables_sva,  "DEGs_SVA_BE_GSE58831.csv"), sep = ",")
deg_gse58831_sva_leek <- read.csv(file.path(paths$tables_sva,  "DEGs_SVA_LEEK_GSE58831.csv"), sep = ",")

deg_gse114922_sva_be <- read.csv(file.path(paths$tables_sva,  "DEGs_SVA_BE_GSE114922.csv"), sep = ",")


# =========================================================
#                  Significant DEGs 
# =========================================================

logcf_cutoff <- 1
padj_cutoff  <- 0.05

# GSE19429
sig_deg_gse19429_be <- subset(deg_gse19429_sva_be, adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff)
sig_deg_gse19429_be <- sig_deg_gse19429_be[!(is.na(sig_deg_gse19429_be$SYMBOL)), ]
sig_deg_gse19429_leek <- subset(deg_gse19429_sva_leek, adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff)
sig_deg_gse19429_leek <- sig_deg_gse19429_leek[!(is.na(sig_deg_gse19429_leek$SYMBOL)), ]

# GSE58831
sig_deg_gse58831_be <- subset(deg_gse58831_sva_be, adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff)
sig_deg_gse58831_be <- sig_deg_gse58831_be[!(is.na(sig_deg_gse58831_be$SYMBOL)), ]
sig_deg_gse58831_leek <- subset(deg_gse58831_sva_leek, adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff)
sig_deg_gse58831_leek <- sig_deg_gse58831_leek[!(is.na(sig_deg_gse58831_leek$SYMBOL)), ]

# GSE114922
sig_deg_gse114922_be <- subset(deg_gse114922_sva_be, padj < padj_cutoff & abs(log2FoldChange) > logcf_cutoff)
sig_deg_gse114922_be <- sig_deg_gse114922_be[!(is.na(sig_deg_gse114922_be$SYMBOL)), ]





# =========================================================
#       Summary of significant DEGs across datasets (SVA)
# =========================================================

deg_summary_df <- data.frame(
  dataset = c(
    rep("GSE19429", 4),
    rep("GSE58831", 4),
    rep("GSE114922", 2)
  ),
  method = c(
    rep("BE", 2), rep("Leek", 2),
    rep("BE", 2), rep("Leek", 2),
    rep("BE", 2)
  ),
  direction = c(
    rep(c("Up in MDS", "Down in MDS"), 5)
  ),
  n_genes = c(
    # GSE19429
    sum(sig_deg_gse19429_be$logFC > 0),
    sum(sig_deg_gse19429_be$logFC < 0),
    sum(sig_deg_gse19429_leek$logFC > 0),
    sum(sig_deg_gse19429_leek$logFC < 0),
    
    # GSE58831
    sum(sig_deg_gse58831_be$logFC > 0),
    sum(sig_deg_gse58831_be$logFC < 0),
    sum(sig_deg_gse58831_leek$logFC > 0),
    sum(sig_deg_gse58831_leek$logFC < 0),
    
    # GSE114922 (Leek not available - convergence failed)
    sum(sig_deg_gse114922_be$log2FoldChange > 0),
    sum(sig_deg_gse114922_be$log2FoldChange < 0)
  )
)

# Explicit ordering for consistent facet and x-axis display
deg_summary_df$dataset <- factor(
  deg_summary_df$dataset,
  levels = c("GSE19429", "GSE58831", "GSE114922")
)
deg_summary_df$method <- factor(
  deg_summary_df$method,
  levels = c("BE", "Leek")
)

deg_summary_plot <- ggplot(deg_summary_df, aes(
  x = method,
  y = n_genes,
  fill = direction
)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = scales::comma(n_genes)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 4
  ) +
  facet_wrap(~ dataset, nrow = 1, scales = "free_x") +
  scale_fill_manual(
    values = c(
      "Up in MDS" = "#CE2915",
      "Down in MDS" = "#0096FF"
    )
  ) +
  scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.15))
  ) +
  labs(
    title = "Significant DEGs across MDS cohorts after SVA correction",
    subtitle = "DEGs defined as adjusted P-value < 0.05 and |log2FC| > 1. Leek method not available for GSE114922 (convergence failure).",
    x = "SVA method",
    y = "Number of significant genes",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10),
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey90", color = NA),
    axis.text.x = element_text(size = 11)
  )

save_plot(deg_summary_plot, filename = "cross_dataset_deg_summary_sva_bar.png", dir = paths$results_sva,
          width = 11, height = 7)



# =========================================================
#                  Save session information
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(
  filename = "sessionInfo_cross_dataset_comparison_sva.txt",
  dir = paths$logs_sva,
  label = "Cross-dataset transcriptomic comparison of SVA-corrected DEGs(GSE19429, GSE58831, GSE114922)"
)

message("[OUTPUT] Session information saved to: ", paths$logs_sva)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Cross-dataset comparison of SVA-corrected DEGs completed.")
message("=================================================")
