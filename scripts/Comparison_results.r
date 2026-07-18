# =========================================================
# Project      : MDS Genomics - Cross-Dataset Comparison
# Datasets     : GSE19429, GSE58831 (microarray) + GSE114922 (RNA-seq)
# Script       : Cross-dataset comparison of Hallmark H (GSEA) and 
#                KEGG (ORA) results
# Description  : Cross-dataset comparison of pathway enrichment and DEG
#                concordance across independent MDS transcriptomic cohorts.
#                Designed for reproducible pre/post SVA correction comparison.
# =========================================================

# =========================================================
#                   Functions
# =========================================================

# 02-ComplexHeatmap function
heatmap_plot <- function(nes_matrix, sig_labels, top_annotation = NULL,
                    cluster_rows = TRUE, cluster_columns = FALSE, title){

  message("[PLOT] Generating Pathways comparison heatmap...")

  # Define symmetric color scale centered on NES = 0

   lim <- max(abs(nes_matrix), na.rm = TRUE)

  col_fun <- circlize::colorRamp2(
    c(-lim, 0, lim),
    c("#4575B4", "white", "#D73027")
  )

  # Generate heatmap

  plot_ht <- ComplexHeatmap::Heatmap(
    nes_matrix,

    name = "NES",

    col = col_fun,

    na_col = "grey90",

    top_annotation = top_annotation,

    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,

    row_names_side = "left",
    row_names_gp = grid::gpar(fontsize = 10),

    column_names_gp = grid::gpar(
      fontsize = 11,
      fontface = "bold"
    ),

    rect_gp = grid::gpar(
      col = "white",
      lwd = 1
    ),

    heatmap_legend_param = list(
      title = "NES",
      title_gp = grid::gpar(fontface = "bold"),
      labels_gp = grid::gpar(fontsize = 10)
    ),

    # Display significance labels (*, **, ***)

    cell_fun = function(j, i, x, y, width, height, fill) {

      grid::grid.text(
        sig_labels[i, j],
        x = x,
        y = y,
        gp = grid::gpar(fontsize = 9)
      )

    },

    column_title = title,
    column_title_gp = grid::gpar(
      fontsize = 14,
      fontface = "bold"
    ),

    row_title = "Pathways",
    row_title_gp = grid::gpar(fontface = "bold")
  )

  message("[PLOT] Heatmap successfully generated.")

  return(plot_ht)

}

# 01-Prepare clean data for plot
prepare_gsea_comparison <- function(gsea_list, min_datasets = 1, min_present  = 3, top_n = NULL, label = "Pathway") {
  
  n_datasets_total <- length(gsea_list)
  
  # Default: require pathway presence in all datasets to generate an NA-free matrix 
  combine_one <- function(df, dataset_label) {
    df <- as.data.frame(df)
    df$dataset <- dataset_label
    df[, c("pathway", "NES", "padj", "dataset")]
  }
  
  message(sprintf("[COMPARE] Preparing %s comparison matrices...", label))
  
  # Combine GSEA results from all datasets into a single long-format table
  pathway_all <- dplyr::bind_rows(
    lapply(names(gsea_list), function(nm)
      combine_one(gsea_list[[nm]], nm))
  )
  
  # Standardize MSigDB Hallmark pathway names for visualization
  if (any(grepl("^HALLMARK_", pathway_all$pathway))) {
    message("[COMPARE] Hallmark pathways detected. Cleaning pathway names...")
    pathway_all$pathway_clean <- tolower(
      gsub("_", " ", sub("^HALLMARK_", "", pathway_all$pathway))
    )
  } else {
    pathway_all$pathway_clean <- pathway_all$pathway
  }
  
  message(sprintf("[COMPARE] Combined %d pathway x dataset rows.", nrow(pathway_all)))
  
  # =======================================================
  # Compute summary statistics for each pathway
  # =======================================================
  
  pathway_stats <- pathway_all %>%
    dplyr::group_by(pathway_clean) %>%
    dplyr::summarise(
      n_present    = sum(!is.na(NES)),
      n_sig        = sum(padj < 0.05, na.rm = TRUE),
      mean_abs_nes = mean(abs(NES), na.rm = TRUE),
      .groups = "drop"
    )
  
  # =======================================================
  # STEP 1 - Filter pathways based on presence across datasets
  # This step controls missing values (NA) in the final matrix
  # =======================================================
  
  present_pool <- pathway_stats[pathway_stats$n_present >= min_present, ]
  
  message(sprintf(
    "[COMPARE] %d / %d %s pathways present in >=%d/%d datasets (NA-free pool).",
    nrow(present_pool), nrow(pathway_stats), label, min_present, n_datasets_total
  ))
  
  # =======================================================
  # STEP 2 - Select pathways for visualization
  # =======================================================
  
  if (!is.null(top_n)) {
    
    # Rank pathways by number of significant datasets, followed by average NES magnitude
    ranking <- present_pool[order(-present_pool$n_sig, -present_pool$mean_abs_nes), ]
    sig_pathways <- head(ranking$pathway_clean, top_n)
    
    n_top_never_sig <- sum(ranking$n_sig[ranking$pathway_clean %in% sig_pathways] == 0)
    if (n_top_never_sig > 0) {
      message(sprintf(
        "[COMPARE] Note: %d of the top %d %s pathways are not significant (padj>=0.05) in any dataset — included only to match the requested panel size.",
        n_top_never_sig, top_n, label
      ))
    }
    
  } else {
    
    
    sig_pathways <- present_pool$pathway_clean[present_pool$n_sig >= min_datasets]
  }
  
  message(sprintf("[COMPARE] %d %s pathways selected for plotting.", length(sig_pathways), label))
  
  # Subset pathways selected for visualization
  plot_df <- pathway_all[pathway_all$pathway_clean %in% sig_pathways, ]
  
  # Build NES matrix
  nes_matrix <- reshape2::dcast(plot_df, pathway_clean ~ dataset, value.var = "NES")
  
  rownames(nes_matrix) <- nes_matrix$pathway_clean
  nes_matrix$pathway_clean <- NULL
  nes_matrix <- as.matrix(nes_matrix)
  
  # Build adjusted p-value matrix
  padj_matrix <- reshape2::dcast(plot_df, pathway_clean ~ dataset, value.var = "padj")
  
  rownames(padj_matrix) <- padj_matrix$pathway_clean
  padj_matrix$pathway_clean <- NULL
  padj_matrix <- as.matrix(padj_matrix)
  
  # Ensure pathway ordering is identical between NES and adjusted p-value matrices
  padj_matrix <- padj_matrix[rownames(nes_matrix), , drop = FALSE]
  
  # Generate significance labels for heatmap annotation
  sig_labels <- matrix("", nrow = nrow(padj_matrix), ncol = ncol(padj_matrix),
                       dimnames = dimnames(padj_matrix))
  sig_labels[padj_matrix < 0.05]  <- "*"
  sig_labels[padj_matrix < 0.01]  <- "**"
  sig_labels[padj_matrix < 0.001] <- "***"
  sig_labels[is.na(padj_matrix)]  <- ""
  
  # Report missing values in the final NES matrix
  n_na <- sum(is.na(nes_matrix))
  message(sprintf("[COMPARE] %d NA cell(s) in final NES matrix.", n_na))
  
  message("[COMPARE] Comparison matrices successfully generated.")
  
  return(list(
    nes_matrix    = nes_matrix,
    padj_matrix   = padj_matrix,
    sig_labels    = sig_labels,
    pathway_all   = pathway_all,
    pathway_stats = pathway_stats
  ))
}




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

message("[LOAD] Reading results from GSEA... (Hallmakrs Pathways - Collection H)")

hallmarks_gse19429  <- read.csv(file.path(paths$tables_gsea,  "results_GSE19429_GSEA_hallmarks_h_df.csv"), sep = ",")
hallmarks_gse58831  <- read.csv(file.path(paths$tables_gsea,  "results_GSE58831_GSEA_hallmarks_h_df.csv"), sep = ",")
hallmarks_gse114922 <- read.csv(file.path(paths$tables_gsea,  "results_GSE114922_GSEA_hallmarks_h_df.csv"), sep = ",")

message("[LOAD] Reading results from GSEA... (KEGG Pathways)")

kegg_gse19429  <- read.csv(file.path(paths$tables_gsea,  "results_GSE19429_GSEA_kegg_df.csv"), sep = ",")
kegg_gse58831  <- read.csv(file.path(paths$tables_gsea,  "results_GSE58831_GSEA_kegg_df.csv"), sep = ",")
kegg_gse114922 <- read.csv(file.path(paths$tables_gsea,  "results_GSE114922_GSEA_kegg_df.csv"), sep = ",")

message("[LOAD] Reading DEG results for 3 different datasets...")

deg_table_gse19429 <- read.csv(file = file.path(paths$tables, "DEGs_GSE19429.csv"), sep = "," )
deg_table_gse58831 <- read.csv(file = file.path(paths$tables, "DEGs_GSE58831.csv"), sep = "," )
deg_table_gse114922 <- read.csv(file = file.path(paths$tables, "DEGs_GSE114922.csv"), sep = "," )


# =========================================================
#        Prepare Data for comparison - Heatmap
# =========================================================

# Hallmarks Pathways - Collection H

res_hallmark <- prepare_gsea_comparison(
  gsea_list = list(
    GSE19429   = hallmarks_gse19429,
    GSE58831   = hallmarks_gse58831,
    GSE114922  = hallmarks_gse114922
  ),
  label = "Hallmark H"
)

# Kegg Pathways

prepare_kegg_gsea <- function(kegg_df) {
  
  message("[PREP] Standardizing KEGG GSEA results...")
  
  kegg_df_clean <- as.data.frame(kegg_df)
  
  colnames(kegg_df_clean)[colnames(kegg_df_clean) == "Description"] <- "pathway"
  colnames(kegg_df_clean)[colnames(kegg_df_clean) == "p.adjust"]    <- "padj"
  
  message("[PREP] Column names successfully standardized.")
  
  return(kegg_df_clean)
}

# Standardize KEGG output column names

kegg_gse19429  <- prepare_kegg_gsea(kegg_gse19429)
kegg_gse58831  <- prepare_kegg_gsea(kegg_gse58831)
kegg_gse114922 <- prepare_kegg_gsea(kegg_gse114922)


## Prepare Data for comparison (KEGG Pathways)

res_kegg <- prepare_gsea_comparison(
  gsea_list = list(
    GSE19429   = kegg_gse19429,
    GSE58831   = kegg_gse58831,
    GSE114922  = kegg_gse114922
  ),
  top_n = 50,    
  label = "KEGG"
)


# =========================================================
#              Complex Heatmap - Pathways 
# =========================================================

technology <- c(
  GSE114922 = "RNA-seq",
  GSE19429 = "Microarray",
  GSE58831 = "Microarray"
)

col = list(
  Technology = c(
    "RNA-seq"     = "#4DBBD5",
    "Microarray"  = "#E64B35"
  )
)

col_annotation <- HeatmapAnnotation(
  Technology = technology[colnames(res_hallmark[["nes_matrix"]])],
  col = col
)

# Hallmark H enrichment heatmap

open_png(filename = "gsea_heatmap_hallmark_H.png", dir = paths$plots_comparison,
         width = 1600, height = 1800)

heatmap_plot(res_hallmark[["nes_matrix"]], sig_labels = res_hallmark[["sig_labels"]], top_annotation = col_annotation, title = "Comparison of MSigDB Hallmark (H) enrichment across MDS datasets")

close_png()

# KEGG enrichment heatmap

open_png(filename = "gsea_heatmap_kegg.png", dir = paths$plots_comparison,
         width = 1600, height = 1800)

heatmap_plot(res_kegg[["nes_matrix"]], sig_labels = res_kegg[["sig_labels"]], top_annotation = col_annotation, title = "Comparison of Kegg enrichment across MDS datasets")

close_png()


# =========================================================
#          .      Significance DEGs 
# =========================================================

# Identify significant DEGs: GSE19429
logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg_gse19429 <- subset(
  deg_table_gse19429,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

sig_deg_gse19429 <- sig_deg_gse19429[!(is.na(sig_deg_gse19429$SYMBOL) & is.na(sig_deg_gse19429$ENTREZID) & is.na(sig_deg_gse19429$GENENAME)), ]

## Identify significant DEGs: GSE58831
logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg_gse58831 <- subset(
  deg_table_gse58831,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

sig_deg_gse58831 <- sig_deg_gse58831[!(is.na(sig_deg_gse58831$SYMBOL) & is.na(sig_deg_gse58831$ENTREZID) & is.na(sig_deg_gse58831$GENENAME)), ]

## Identify significant DEGs: GSE114922
logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg_gse114922 <- subset(
  deg_table_gse114922,
  padj < padj_cutoff & abs(log2FoldChange) > logcf_cutoff
)

sig_deg_gse114922 <- sig_deg_gse114922[!(is.na(sig_deg_gse114922$SYMBOL) & is.na(sig_deg_gse114922$ENTREZID) & is.na(sig_deg_gse114922$GENENAME)), ]


# ===================================================================
# DEG overlap analysis at gene-symbol level (up/down regulated genes)
# ===================================================================

deg_symbols_up <- list(
  GSE19429   = sig_deg_gse19429$SYMBOL[sig_deg_gse19429$logFC > 0],
  GSE58831   = sig_deg_gse58831$SYMBOL[sig_deg_gse58831$logFC > 0],
  GSE114922  = sig_deg_gse114922$SYMBOL[sig_deg_gse114922$log2FoldChange > 0]
)

deg_symbols_down <- list(
  GSE19429   = sig_deg_gse19429$SYMBOL[sig_deg_gse19429$logFC < 0],
  GSE58831   = sig_deg_gse58831$SYMBOL[sig_deg_gse58831$logFC < 0],
  GSE114922  = sig_deg_gse114922$SYMBOL[sig_deg_gse114922$log2FoldChange < 0]
)

# Upregulated DEG overlap
venn_up <-ggVennDiagram::ggVennDiagram(deg_symbols_up, label = "count", label_alpha = 0) +
  scale_fill_gradient(
    low = "white",
    high = "#C44E52",
    guide = "none"
  ) +
  labs(
    title = "Significant DEGs (up in MDS) - Gene overlap across datasets"
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

# Downregulated DEG overlap
venn_down <-ggVennDiagram::ggVennDiagram(deg_symbols_down, label = "count", label_alpha = 0) +
  scale_fill_gradient(
    low = "white",
    high = "#4C72B0",
    guide = "none"
  ) +
  labs(
    title = "Significant DEGs (down in MDS) - Gene overlap across datasets"
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    )
  )

# Save plots
save_plot(venn_up, filename = "cross_dataset_deg_venn_up.pdf", dir = paths$plots_comparison, bg = "transparent", width = 8, height = 7)
save_plot(venn_down, filename = "cross_dataset_deg_venn_down.pdf", dir = paths$plots_comparison, width = 8, height = 7)

# Lista esplicita dei geni condivisi da tutti e tre - utile da citare nel testo del report
shared_up_all3 <- Reduce(intersect, deg_symbols_up)
shared_down_all3 <- Reduce(intersect, deg_symbols_down)

message(sprintf("[COMPARE] Shared upregulated genes across all datasets: %d -> %s",
                length(shared_up_all3), paste(shared_up_all3, collapse = ", ")))
message(sprintf("[COMPARE] Shared downregulated genes across all datasets: %d -> %s",
                length(shared_down_all3), paste(shared_down_all3, collapse = ", ")))


# =============================================================
#      # Pairwise log2FC concordance analysis between datasets
# =============================================================

#--------------------------------------------------
# Function: calculate log2FC correlation and generate scatter plot
#--------------------------------------------------
plot_logfc_concordance <- function(deg_table_x, x_name, logfc_col_x,
                                   deg_table_y, y_name, logfc_col_y,
                                   padj_col_x, padj_col_y,
                                   plot_name) {

  stopifnot(
    anyDuplicated(deg_table_x$SYMBOL) == 0,
    anyDuplicated(deg_table_y$SYMBOL) == 0
  )
  
  merged <- dplyr::inner_join(
    deg_table_x[, c("SYMBOL", logfc_col_x, padj_col_x)],
    deg_table_y[, c("SYMBOL", logfc_col_y, padj_col_y)],
    by = "SYMBOL"
  )
  
  message("Ok merged.")
  colnames(merged) <- c("SYMBOL", "logFC_x", "padj_x", "logFC_y", "padj_y")
  merged <- merged[!is.na(merged$logFC_x) & !is.na(merged$logFC_y), ]
  
  merged$both_sig <- merged$padj_x < 0.05 & merged$padj_y < 0.05
  
  cor_val <- cor(merged$logFC_x, merged$logFC_y, method = "pearson", use = "complete.obs")
  
  message(sprintf("[COMPARE] %s vs %s: Pearson r = %.3f | %d shared genes (%d significant in both).",
                  x_name, y_name, cor_val, nrow(merged), sum(merged$both_sig)))

  p <- ggplot2::ggplot(merged, ggplot2::aes(x = logFC_x, y = logFC_y, color = both_sig)) +
    ggplot2::geom_point(alpha = 0.4, size = 1) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      color = "grey60"
    ) +
    ggplot2::geom_smooth(
      method = "lm",
      color = "black",
      se = FALSE,
      linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = c(`FALSE` = "grey70", `TRUE` = "#CE2915")) +
    ggplot2::labs(
      title = plot_name,
      x = paste0("log2FC ", x_name),
      y = paste0("log2FC ", y_name),
      color = NULL
    ) +
    ggplot2::theme_classic(base_size = 14)
  
  
  return(list(
    plot = p,
    correlation = cor_val,
    merged = merged
  ))
  
}



#--------------------------------------------------
# Prepare input data: collapse to gene-level
#--------------------------------------------------
# Collapse probe-level results to gene-level representation
collapse_deg <- function(df, padj_col) {
  
  df %>%
    filter(!is.na(SYMBOL)) %>%
    arrange(SYMBOL, .data[[padj_col]]) %>%
    group_by(SYMBOL) %>%
    slice_head(n = 1) %>%
    ungroup()
  
}

deg_table_gse114922_gene <- collapse_deg(deg_table_gse114922, "padj")
deg_table_gse19429_gene <- collapse_deg(deg_table_gse19429, "adj.P.Val")
deg_table_gse58831_gene <- collapse_deg(deg_table_gse58831, "adj.P.Val")

#--------------------------------------------------
# Cross-cohort concordance: Affymetrix arrays
#--------------------------------------------------

cor_array_array <- plot_logfc_concordance(
  deg_table_gse19429_gene, "GSE19429", "logFC",
  deg_table_gse58831_gene, "GSE58831", "logFC",
  "adj.P.Val", "adj.P.Val",
  plot_name =  "Cross-cohort log2FC concordance (Affymetrix arrays)"
)

#--------------------------------------------------
# Cross-platform concordance: RNA-seq vs Affymetrix array
#--------------------------------------------------

cor_rnaseq_array1 <- plot_logfc_concordance(
  deg_table_gse19429_gene, "GSE19429", "logFC",
  deg_table_gse114922_gene, "GSE114922", "log2FoldChange",
  "adj.P.Val", "padj",
  plot_name = "Cross-platform log2FC concordance (GSE114922 RNA-seq vs GSE19429 Affymetrix)"
)

#--------------------------------------------------
# Cross-platform concordance: RNA-seq vs Affymetrix array
#--------------------------------------------------

cor_rnaseq_array2 <- plot_logfc_concordance(
  deg_table_gse58831_gene, "GSE58831", "logFC",
  deg_table_gse114922_gene, "GSE114922", "log2FoldChange",
  "adj.P.Val", "padj",
  plot_name = "Cross-platform log2FC concordance (GSE114922 RNA-seq vs GSE58831 Affymetrix)"
)

# Saving plots

save_plot(plot = cor_array_array[["plot"]], filename = "log2FC_concordance_GSE19429_vs_GSE58831.png", dir = paths$plots_comparison,
          width = 11, height = 9)

save_plot(plot = cor_rnaseq_array1[["plot"]], filename = "log2FC_concordance_GSE114922_vs_GSE19429.png", dir = paths$plots_comparison,
          width = 11, height = 9)
save_plot(plot = cor_rnaseq_array2[["plot"]], filename = "log2FC_concordance_GSE114922_vs_GSE58831.png", dir = paths$plots_comparison,
          width = 11, height = 9)


# =========================================================
#       Summary of significant DEGs across datasets
# =========================================================

deg_summary_df <- data.frame(
  dataset = rep(c("GSE19429", "GSE58831", "GSE114922"), each = 2),
  direction = rep(c("Up in MDS", "Down in MDS"), times = 3),
  n_genes = c(
    sum(sig_deg_gse19429$logFC > 0),  sum(sig_deg_gse19429$logFC < 0),
    sum(sig_deg_gse58831$logFC > 0),  sum(sig_deg_gse58831$logFC < 0),
    sum(sig_deg_gse114922$log2FoldChange > 0), sum(sig_deg_gse114922$log2FoldChange < 0)
  )
)

deg_summary_plot <- ggplot(deg_summary_df,aes(
    x = dataset,
    y = n_genes,
    fill = direction
  )
) +
  geom_col(
    position = ggplot2::position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = scales::comma(n_genes)),
    position = position_dodge(width = 0.8),
    vjust = -0.4,
    size = 4
  ) +
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
    title = "Significant differentially expressed genes across MDS cohorts",
    subtitle = "DEGs defined as adjusted P-value < 0.05 and |log2FC| > 1",
    x = NULL,
    y = "Number of significant genes",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 11
    ),
    legend.position = "top",
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

save_plot(deg_summary_plot, filename = "cross_dataset_deg_summary_bar.png", dir = paths$plots_comparison,
          width = 9, height = 8)


# =========================================================
#                  Save session information
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(
  filename = "sessionInfo_cross_dataset_comparison.txt",
  dir = paths$logs,
  label = "Cross-dataset transcriptomic comparison (GSE19429, GSE58831, GSE114922)"
)

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Cross-dataset comparison completed successfully.")
message("[PIPELINE] Generated pathway enrichment heatmaps, DEG overlap analysis, and log2FC concordance plots.")
message("=================================================")
