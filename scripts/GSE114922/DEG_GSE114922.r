# =========================================================
# Project      : MDS Genomics - Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE114922 (Illumina HiSeq 4000)
# Samples      : 82 MDS patients + 8 healthy controls, CD34+ bone marrow HSC
# Script       : Quality Control, Exploratory Data Analysis and Differential Expression Analysis
# Description  : End-to-end workflow: Build DESeq2 object and QC
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("SetupEnvironment/00_paths.R")
source("SetupEnvironment/01_environment.R")
source("SetupEnvironment/02_seed.R")
source("SetupEnvironment/03_helper_functions.R")

set_seed(1234)

message("=== Starting GSE114922 RNA-seq preprocessing and QC ===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Loading raw count matrix...")

counts_raw <- read.delim(gzfile("GSE114922_raw_counts_GRCh38.p13_NCBI.tsv.gz"), row.names = 1)

# Download GEO metadata

gse_id <- "GSE114922"
gse_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

if (is.list(gse_list) && !is(gse_list, "ExpressionSet")) {
  message("getGEO() returned ", length(gse_list), " ExpressionSet object(s).")
  gse <- gse_list[[1]]
} else {
  gse <- gse_list
}

meta <- pData(gse)
rownames(meta) <- meta$geo_accession


# =========================================================
#                  Metadata cleaning 
# =========================================================

status_col <- grep("disease status", colnames(meta), ignore.case = TRUE, value = TRUE)

if (length(status_col) == 0) {
  stop(
    "Could not find a 'disease status' column in pData(gse). ",
    "Inspect colnames(pData(gse)) manually and update this script:\n",
    paste(colnames(meta), collapse = "\n")
  )
} else message("[META] Disease status column identified.")

if (length(status_col) > 1) {
  message("Multiple candidate columns found, using the first: ", status_col[1])
}
status_col <- status_col[1]

group_raw <- as.character(meta[[status_col]])
message("Raw values found in '", status_col, "':")
print(table(group_raw, useNA = "ifany"))

# Recode column metadata into a clean two-level factor

group <- rep(NA_character_, length(group_raw))
group[grepl("MDS|myelodysplastic", group_raw, ignore.case = TRUE)] <- "MDS"
group[grepl("healthy|control", group_raw, ignore.case = TRUE)] <- "HC"
group <- factor(group, levels = c("HC", "MDS"))

if (any(is.na(group))) {
  stop("Some samples could not be classified into MDS/HC. Check group_raw values.")
}

pData(gse)$group <- group
message("[META] Sample distribution:")
print(table(pData(gse)$group))

names(group) <- meta$geo_accession

meta$group <- group


# =========================================================
#                 subset only CD34 + cells  
# =========================================================

cd34_idx <- grepl("CD34\\+ hematopoietic", meta$characteristics_ch1.1)

meta_cd34 <- meta[cd34_idx, ]
counts_cd34 <- counts_raw[, meta_cd34$geo_accession, drop = FALSE]
meta_cd34$group <- group[rownames(meta_cd34)]

# Consistency checks

stopifnot(all(ncol(counts_cd34) == nrow(meta_cd34)))
stopifnot(ncol(counts_cd34) > 1)
stopifnot(nrow(counts_cd34) > 1000)

stopifnot(!any(duplicated(rownames(counts_cd34))))
stopifnot(!any(duplicated(colnames(counts_cd34))))

stopifnot(!any(is.na(counts_cd34)))
stopifnot(all(is.finite(as.matrix(counts_cd34))))

message("[LOAD] Count matrix dimensions: ", nrow(counts_cd34), " genes x ", ncol(counts_cd34), " samples")

# Verify sample order

stopifnot(all(colnames(counts_cd34) == meta_cd34$geo_accession))

# Build DESeq2 sample metadata

colData <- data.frame(
  row.names = meta_cd34$geo_accession,
  condition = factor(meta_cd34$group, levels = c("HC", "MDS"))
)

message("[META] Sample distribution:")

print(table(colData$condition))

stopifnot(all(colnames(counts_cd34) == rownames(colData)))
stopifnot(length(unique(colData$condition)) >= 2)

stopifnot(all(colnames(counts_cd34) %in% rownames(colData)))
stopifnot(all(colnames(counts_cd34) == rownames(colData)))


# =========================================================
#                   DESeq2 
# =========================================================

# Build a DESeq Dataset object

dds <- DESeqDataSetFromMatrix(
  countData = counts_cd34,
  colData = colData,
  design = ~ condition
  )

#  Filter lowly expressed genes

dds <- dds[rowSums(counts(dds)) > 10, ]

# Set reference level and run DESeq2 analysis

dds$condition <- relevel(dds$condition, ref = "HC")
dds <- DESeq(dds)
res <- results(dds)

message("[DESeq2] Significant DE genes (padj < 0.05): ",
        sum(res$padj < 0.05, na.rm = TRUE))

stopifnot(!is.null(dispersions(dds)))

# Variance stabilization using vst()

vsdata <- vst(dds, blind = FALSE)
vsdata_df <- assay(vsdata)


# =========================================================
#                  Quality Control 
# =========================================================

# Distribution of VST-transformed expression values

group_palette <- setNames(brewer.pal(3, "Set1")[1:2], levels(meta_cd34$group))
sample_colors <- group_palette[as.character(meta_cd34$group)]

open_png(filename = "01_boxplot_intensity_distributions.png", dir = paths$plots_gse114922_qc,
         width = 1600, height = 900)
par(mar = c(8, 4, 4, 2))
boxplot(vsdata_df,
        main = "Distribution of VST-transformed expression values per sample (GSE114922)",
        ylab = "VST expression",
        col = sample_colors,
        las = 2, cex.axis = 0.5, outline = FALSE)
legend("topright", legend = levels(meta_cd34$group), fill = group_palette, bty = "n")
close_png()

open_png(filename = "02_density_plot_vst_expression.png", dir = paths$plots_gse114922_qc, 
         width = 1200, height = 800)
plot(density(vsdata_df[, 1]), col = sample_colors[1], lwd = 1,
     main = "Density distribution of VST-transformed expression values across samples",
     xlab = "VST expression", ylim = c(0, 0.8))
for (i in 2:ncol(vsdata_df)) {
  lines(density(vsdata_df[, i]), col = sample_colors[i], lwd = 0.8)
}
legend("topright", legend = levels(meta_cd34$group), col = group_palette, lwd = 2, bty = "n")
close_png()

# Library size distribution

lib_sizes_dds <- colSums(counts(dds))
lib_dds_df <- data.frame(sample = names(lib_sizes_dds), lib_size = lib_sizes_dds, group = meta_cd34$group)

libsize_plot <- ggplot(lib_dds_df, aes(x = reorder(sample, lib_sizes_dds), y = lib_size, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("#CE2915", "#0096FF")) +
  labs(title = "Library size per sample (GSE114922)", x = "", y = "Total mapped reads") +
  theme_bw()

save_plot(libsize_plot, filename = "03_library_size.png", dir = paths$plots_gse114922_qc,
          width = 9, height = 9)


# =========================================================
#               Outlier Detection (QC)
# =========================================================

message("[QC] Starting outlier detection...")

# ---------------------------------------------------------
# Safety checks
# ---------------------------------------------------------

stopifnot(exists("vsdata_df"))
stopifnot(exists("counts_cd34"))
stopifnot(is.matrix(vsdata_df) || is.data.frame(vsdata_df))

vsdata_df <- as.matrix(vsdata_df)

# ---------------------------------------------------------
# Sample-wise median VST expression (robust central tendency)
# ---------------------------------------------------------

message("[QC] Computing sample-wise median VST expression...")

sample_medians <- matrixStats::colMedians(vsdata_df, na.rm = TRUE)

med_center <- median(sample_medians, na.rm = TRUE)

med_iqr <- IQR(sample_medians, na.rm = TRUE)

outliers_median <- names(sample_medians)[abs(sample_medians - med_center) > 3 * med_iqr]

message(sprintf(
  "[QC] Median-VST rule flagged %d sample(s).",
  length(outliers_median)
))

# ---------------------------------------------------------
# Sample-sample correlation QC
# ---------------------------------------------------------

message("[QC] Computing sample-sample Pearson correlation (VST)...")

cor_mat <- cor(vsdata_df, method = "pearson", use = "pairwise.complete.obs")
diag(cor_mat) <- NA

avg_cor <- rowMeans(cor_mat, na.rm = TRUE)

cor_mean <- mean(avg_cor, na.rm = TRUE)
cor_sd <- sd(avg_cor, na.rm = TRUE)

cor_threshold <- cor_mean - 3 * cor_sd

outliers_corr <- names(avg_cor)[avg_cor < cor_threshold]

message(sprintf(
  "[QC] Correlation rule flagged %d sample(s).",
  length(outliers_corr)
))

# ---------------------
# Library size QC 
# ---------------------

message("[QC] Computing library size per sample (raw counts)...")

lib_sizes <- colSums(counts_cd34)

lib_center <- median(lib_sizes, na.rm = TRUE)
lib_iqr <- IQR(lib_sizes, na.rm = TRUE)

outliers_libsize <- names(lib_sizes)[abs(lib_sizes - lib_center) > 3 * lib_iqr]

message(sprintf(
  "[QC] Library-size rule flagged %d sample(s).",
  length(outliers_libsize)
))

# ---------------------------------------------------------
# Combine outliers
# ---------------------------------------------------------

outlier_samples <- union(union(outliers_median, outliers_corr), outliers_libsize)

message(sprintf(
  "[QC] Total unique outliers: %d sample(s).",
  length(outlier_samples)
))

# ---------------------------------------------------------
# Output 
# ---------------------------------------------------------

qc_outliers_df <- data.frame(
  sample = outlier_samples,
  reason_median_vst = outlier_samples %in% outliers_median,
  reason_low_correlation = outlier_samples %in% outliers_corr,
  reason_library_size = outlier_samples %in% outliers_libsize,
  stringsAsFactors = FALSE
)

message("[QC] Outlier detection completed (Not removed).")


# =========================================================
#           Exploratory Data Analysis (EDA) 
# =========================================================

# --------------------------
# PCA (VST transformation)
# --------------------------

pca_res <- prcomp(t(vsdata_df), center = TRUE, scale. = FALSE)

pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = meta_cd34$group,
  outlier = colnames(vsdata_df) %in% outlier_samples
)

pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

# PCA coloured by disease group

pca_plot_group <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = outlier)) +
  geom_point(size = 2.2, alpha = 0.85) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  scale_color_manual(values = c("#CE2915", "#0096FF")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  ggtitle("PCA of GSE114922 samples") +
  theme_bw()

save_plot(pca_plot_group, filename = "03_pca_plot.png", dir = paths$plots_gse114922_eda)

# --------------------------
# Scree plot
# --------------------------

var_df <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  variance = pca_var[1:10] * 100,
  cumvar = cumsum(pca_var[1:10]) * 100
)

open_png(filename = "04_pca_scree_plot.png", dir = paths$plots_gse114922_eda,
         width = 1000, height = 700)

par(mar = c(5, 5, 4, 2))

# Bar plot: percentage of variance explained by each principal component
bar_centers <- barplot(
  var_df$variance,
  names.arg = var_df$PC,
  ylim = c(0, max(var_df$variance) * 1.2),
  ylab = "% variance explained",
  main = "PCA Scree Plot (Top 10 PCs)",
  col = "grey70"
)

par(new = TRUE)

# Overlay cumulative explained variance on a secondary y-axis
plot(
  bar_centers, var_df$cumvar,
  type = "b",
  axes = FALSE,
  xlab = "", ylab = "",
  xlim = par("usr")[1:2],
  ylim = c(0, 100),
  col = "red", pch = 16
)

# Add secondary axis for cumulative explained variance
axis(side = 4)
mtext("Cumulative variance (%)", side = 4, line = 3)

close_png()

# ---------------------------------
# Hierarchical clustering (samples)
# ---------------------------------

message("[QC] Computing sample distance matrix...")

sample_dist <- dist(t(vsdata_df), method = "euclidean")

hc <- hclust(sample_dist, method = "average")

# plot

open_png(filename = "05_hierarchical_clustering.png", dir = paths$plots_gse114922_eda,
         width = 2200, height = 900)

par(mar = c(5, 4, 4, 2))

plot(hc, labels = FALSE,
     main = "Hierarchical clustering of samples (Euclidean, average linkage)",
     xlab = "",
     sub = ""
)

# Color annotation aligned to dendrogram order

ordered_groups <- meta_cd34$group[hc$order]

# Preserve group-to-colour mapping after dendrogram reordering

ordered_colors <- group_palette[as.character(ordered_groups)]

rect_positions <- seq_along(hc$order)

# Add group annotation below the dendrogram

for (i in seq_along(hc$order)) {
  axis(1, at = i, labels = FALSE, col.ticks = ordered_colors[i], lwd.ticks = 3)
}

# legend
legend("top", legend = names(group_palette), fill = group_palette, bty = "n")

close_png()

message("[QC] Hierarchical clustering plot saved.")


# ---------------------------------
# Sample-sample correlation heatmap
# ---------------------------------

message("[QC] Computing sample-sample correlation matrix...")

sample_cor <- cor(vsdata_df, method = "pearson", use = "pairwise.complete.obs")

# Sample annotations

annotation_col <- data.frame(Group = meta_cd34$group)
rownames(annotation_col) <- colnames(vsdata_df)

# Ensure annotation order matches correlation matrix

annotation_col <- annotation_col[colnames(sample_cor), , drop = FALSE]

stopifnot(identical(rownames(annotation_col), colnames(cor_mat)))
stopifnot(all(rownames(annotation_col) == colnames(cor_mat)))

open_png(filename = "06_sample_correlation_heatmap.png", dir = paths$plots_gse114922_eda,
         width = 1400, height = 1400)

pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_palette),
  show_rownames = FALSE, 
  show_colnames = FALSE,
  main = "Sample-sample Pearson correlation (GSE114922)"
)

close_png()


# ---------------------------------
# Sample-distance heatmap (1 - correlation)
# ---------------------------------

message("[QC] Computing sample distance matrix (1 - cor)...")

sample_dist_1mcor <- as.dist(1 - sample_cor)

open_png(filename = "07_sample_distance_heatmap.png", dir = paths$plots_gse114922_eda,
         width = 1400, height = 1400)
pheatmap(
  as.matrix(sample_dist_1mcor),
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_palette),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Sample distance heatmap (1 - Pearson correlation)"
)

close_png()


# =========================================================
#           DIFFERENTIAL EXPRESSION ANALYSIS (DEGs)
# =========================================================

deg_table_raw <- results(dds, contrast = c("condition", "MDS", "HC"), alpha = 0.05)
deg_table <- as.data.frame(deg_table_raw)
deg_table$ENTREZID <- rownames(deg_table)

deg_table$SYMBOL <- mapIds(
  org.Hs.eg.db,
  keys = rownames(deg_table),
  keytype = "ENTREZID",
  column = "SYMBOL",
  multiVals = "first"
)

deg_table$GENENAME <- mapIds(
  org.Hs.eg.db,
  keys = rownames(deg_table),
  keytype = "ENTREZID",
  column = "GENENAME",
  multiVals = "first"
)

# ------------------------------------------
# Significant DEG summary
# ------------------------------------------

logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg <- subset(
  deg_table,
  padj < padj_cutoff & abs(log2FoldChange) > logcf_cutoff
)

sig_deg <- sig_deg[!is.na(sig_deg$padj),]
n_before <- nrow(sig_deg)

sig_deg <- sig_deg[!(is.na(sig_deg$SYMBOL) & is.na(sig_deg$ENTREZID) & is.na(sig_deg$GENENAME)), ]

message(sprintf(
  "[DEG] Removed %d significant gene(s) lacking all annotation fields (SYMBOL, ENTREZID and GENENAME).",
  n_before - nrow(sig_deg)
))

message(sprintf(
  "[DEG] Significant genes (FDR < %.2f, |log2FC| > %.1f): %d upregulated, %d downregulated (%d total).",
  padj_cutoff,
  logcf_cutoff,
  sum(sig_deg$log2FoldChange > 0),
  sum(sig_deg$log2FoldChange < 0),
  nrow(sig_deg)
))


# =========================================================
#           Plotting DEGs
# =========================================================

# ------------------------------------------
# Volcano Plot
# ------------------------------------------

message("[DEG] Generating volcano plot...")

deg_table$significance <- "Not Significant"
deg_table$significance[deg_table$padj < padj_cutoff & deg_table$log2FoldChange > logcf_cutoff] <- "Up in MDS"
deg_table$significance[deg_table$padj < padj_cutoff & deg_table$log2FoldChange < -logcf_cutoff] <- "Down in MDS"

volcano_plot_deg <- ggplot(deg_table,aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_vline(
    xintercept = c(-logcf_cutoff, logcf_cutoff),
    linetype = "dashed",
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    color = "grey40"
  ) +
  scale_color_manual(
    values = c(
      "Up in MDS" = "#CE2915",
      "Down in MDS" = "#0096FF",
      "Not Significant" = "grey75"
    )
  ) +
  labs(
    title = "Volcano plot of differential expression (GSE114922)",
    x = expression(log[2] ~ fold~change),
    y = expression(-log[10] ~ adjusted~italic(P)),
    color = NULL
  ) +
  theme_bw()

save_plot(volcano_plot_deg, filename = "08_volcano_plot.png", dir = paths$plots_gse114922_deg, 
          width = 8, height = 6)

# ------------------------------------------
# MA Plot
# ------------------------------------------

message("[DEG] Generating MA plot...")

open_png(filename = "09_MA_plot_MDS_vs_HC.png", dir = paths$plots_gse114922_deg,
         width = 1000, height = 800)

plotMA(
  deg_table_raw,
  ylim = c(-7, 7),
  alpha = 0.05,
  main = "MA plot: MDS vs Healthy Controls"
)
abline(h = c(-1, 1), col = "grey50", lty = 2)

close_png()

# -------------------
# Heatmap
# ------------------

# Select top DEGs

top_degs <- 50

# Filter significant DEGs and order by padj

deg_annotated <- sig_deg[!is.na(sig_deg$SYMBOL) & sig_deg$SYMBOL!="", ]
deg_annotated <- deg_annotated[order(deg_annotated$padj), ]
deg_annotated <- deg_annotated[!duplicated(deg_annotated$SYMBOL), ]

message(sprintf(
  "[DEG] %d significant genes available for heatmap.",
  nrow(deg_annotated)
))

# Sanity check

if (nrow(deg_annotated) < top_degs) {
  
  warning(sprintf(
    "[DEG] Only %d significant genes available (requested %d). Heatmap will show all genes.",
    nrow(deg_annotated),
    top_degs
  ))
  
  top_degs <- nrow(deg_annotated)
}

# Select top genes

top_genes <- head(rownames(deg_annotated), top_degs)

# Extract VST expression matrix

heatmap_matrix <- vsdata_df[top_genes, , drop = FALSE]

# Row-wise Z-score normalization

heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

# Row annotation

if (!is.null(sig_deg$SYMBOL)) {
  
  gene_labels <- sig_deg$SYMBOL[
    match(
      rownames(heatmap_matrix_scaled),
      rownames(sig_deg)
    )
  ]
  
  rownames(heatmap_matrix_scaled) <- make.unique(gene_labels)
  
}

# Sample annotation

ann_col <- HeatmapAnnotation(
  Group = colData$condition,
  col = list(Group = group_palette)
)

# Plot
open_png(filename = "10-Heatmap_TopDegs_MDS_vs_HC.png", dir = paths$plots_gse114922_deg,
         width = 1600, height = 1800)

ComplexHeatmap::Heatmap(
  heatmap_matrix_scaled,
  name = "Z-score",
  top_annotation = ann_col,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_title = "Top 50 Differentially Expressed Genes (MDS vs Healthy Controls)",
  row_title = sprintf("Top %d DEGs", top_degs),
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Row\nZ-score"
  )
)

close_png()


# =========================================================
#                Save outputs
# =========================================================

message("[OUTPUT] Saving DESeq2 objects...")

save_rds(dds, filename = "DESeq2_object.rds", dir = paths$results)

message("[OUTPUT] Saving differential expression results...")

save_csv(deg_table, filename = "DEGs_GSE114922.csv", dir = paths$tables)

message("[OUTPUT] Saving background gene list for enrichment analyses...")

background <- rownames(dds)

save_csv(background, filename = "background_GSE114922.csv", dir = paths$tables)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_DESeq2_object.txt", dir = paths$logs, label = "Build DESeq2 object - GSE114922")

message("[OUTPUT] Session information saved to:" , paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Differential expression analysis completed successfully for GSE114922 (MDS vs Healthy Controls).")
message("=================================================")

