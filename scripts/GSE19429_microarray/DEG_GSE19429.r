# =========================================================
# Project      : MDS Genomics - Microarray Analysis Pipeline
# Dataset      : GSE19429 (Affymetrix HG-U133 Plus 2.0, GPL570)
# Samples      : 183 MDS patients + 17 healthy controls, CD34+ bone marrow HSC
# Script       : Quality Control, Exploratory Data Analysis and Differential Expression Analysis
# Description  : End-to-end workflow: QC -> EDA -> Differential Expression
#                Analysis (DEG), comparison MDS vs Healthy Controls.
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("SetupEnvironment/00_paths.R")
source("SetupEnvironment/01_environment.R")
source("SetupEnvironment/02_seed.R")
source("SetupEnvironment/03_helper_functions.R")

set_seed(1234)

message("=== Starting GSE19429 microarray analysis pipeline ===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Downloading normalized expression data from GEO...")

gse_id <- "GSE19429"

gse_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

# check the length of GSE list 

if (is.list(gse_list) && !is(gse_list, "ExpressionSet")) {
  message("getGEO() returned ", length(gse_list), " ExpressionSet object(s).")
  gse <- gse_list[[1]]
} else {
  gse <- gse_list
}

stopifnot(is(gse, "ExpressionSet"))
message("Loaded ExpressionSet: ", ncol(gse), " samples x ", nrow(gse), " probes")


# =========================================================
#                  Metadata cleaning 
# =========================================================

meta <- pData(gse)
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

# Disease subtype (FAB subtype) - used for EDA color coding

subtype_col <- grep("disease subtype", colnames(meta), ignore.case = TRUE, value = TRUE)

if (length(subtype_col) == 1) {
  
  subtype_raw <- as.character(meta[[subtype_col]])
  subtype_clean <- sub("^disease subtype:\\s*", "", subtype_raw)
  subtype_clean <- ifelse(is.na(subtype_clean) | subtype_clean == "", "HC", subtype_clean)
  subtype_clean <- factor(subtype_clean)
  pData(gse)$subtype <- factor(subtype_clean)
  
  message("[META] Disease subtype column processed: ", subtype_col)
  print(table(pData(gse)$subtype, useNA = "ifany"))
  
} else {
  
  message("No disease subtype column found - subtype-based EDA will be skipped.")
  pData(gse)$subtype <- factor("Unknown", levels = "Unknown")
  
}

subtype <- pData(gse)$subtype

# =========================================================
#                  Quality Control 
# =========================================================

expr_raw <- exprs(gse)

message("[QC] GEO preprocessing: ", unique(meta$data_processing))

# Check expression scale and apply log2 transformation if required

max_val_bt <- max(expr_raw, na.rm = TRUE)
min_val_bt <- min(expr_raw, na.rm = TRUE)
message(sprintf("[QC] Expression range before transformation: [%.3f, %.3f]", min_val_bt, max_val_bt))

if (max_val_bt > 50) {
  
  message("Values appear to be on a linear scale -> applying log2(x + 1) transform.")

    if (min_val_bt <= 0) {
      
    warning("Non-positive intensity values detected (", sum(expr_raw <= 0),
            " entries) - flooring to 1e-3 before log2 transform.")
      
    expr_raw[expr_raw <= 0] <- 1e-3
    }
  
  expr_log <- log2(expr_raw + 1)
  
} else {
  message("Values appear to already be log2-scale -> no transform applied.")
  expr_log <- expr_raw
}

max_val_at <- max(expr_log, na.rm = TRUE)
min_val_at <- min(expr_log, na.rm = TRUE)
message(sprintf("[QC] Expression range after transformation: [%.3f, %.3f]", min_val_at, max_val_at))

# Remove Affymetrix control probes

affx_probe <- grepl("^AFFX", rownames(expr_log))
message(sprintf("Removing %d AFFX control probes out of %d total.",sum(affx_probe), nrow(expr_log)))

expr_log <- expr_log[!affx_probe, ]

# Box plot and density plot

group_palette <- setNames(brewer.pal(3, "Set1")[1:2], levels(group))
sample_colors <- group_palette[as.character(group)]

open_png(filename = "01_boxplot_intensity_distributions.png", dir = paths$plots_gse19429_qc,
         width = 1600, height = 900)
par(mar = c(8, 4, 4, 2))
boxplot(expr_log,
        main = "Log2 intensity distribution per sample (GSE19429)",
        ylab = "log2 intensity",
        col = sample_colors,
        las = 2, cex.axis = 0.5, outline = FALSE)
legend("topright", legend = levels(group), fill = group_palette, bty = "n")
close_png()

open_png(filename = "02_density_plot_intensity.png", dir = paths$plots_gse19429_qc, 
         width = 1200, height = 800)
plot(density(expr_log[, 1]), col = sample_colors[1], lwd = 1,
     main = "Density of log2 intensities across samples",
     xlab = "log2 intensity", ylim = c(0, 0.35))
for (i in 2:ncol(expr_log)) {
  lines(density(expr_log[, i]), col = sample_colors[i], lwd = 0.8)
}
legend("topright", legend = levels(group), col = group_palette, lwd = 2, bty = "n")
close_png()


# =========================================================
#               Outlier Detection (QC)
# =========================================================

message("[QC] Starting outlier detection...")

# ---------------------------------------------------------
# 0. Safety checks
# ---------------------------------------------------------

stopifnot(exists("expr_log"))
stopifnot(is.matrix(expr_log) || is.data.frame(expr_log))

expr_log <- as.matrix(expr_log)

# ---------------------------------------------------------
# Sample-wise median intensity (robust central tendency)
# ---------------------------------------------------------

message("[QC] Computing sample-wise median intensities...")

sample_medians <- matrixStats::colMedians(expr_log, na.rm = TRUE)

med_center <- median(sample_medians, na.rm = TRUE)

med_iqr <- IQR(sample_medians, na.rm = TRUE)

outliers_median <- names(sample_medians)[
  abs(sample_medians - med_center) > 3 * med_iqr
]

message(sprintf(
  "[QC] Median-intensity rule flagged %d sample(s).",
  length(outliers_median)
))

# ---------------------------------------------------------
# Sample-sample correlation QC
# ---------------------------------------------------------

message("[QC] Computing sample-sample Pearson correlation...")

cor_mat <- cor(expr_log, method = "pearson", use = "pairwise.complete.obs")
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

# ---------------------------------------------------------
# Combine outliers
# ---------------------------------------------------------

outlier_samples <- union(outliers_median, outliers_corr)

message(sprintf(
  "[QC] Total unique outliers: %d sample(s).",
  length(outlier_samples)
))

# ---------------------------------------------------------
# Output 
# ---------------------------------------------------------

qc_outliers_df <- data.frame(
  sample = outlier_samples,
  reason_median_intensity = outlier_samples %in% outliers_median,
  reason_low_correlation = outlier_samples %in% outliers_corr,
  stringsAsFactors = FALSE
)

message("[QC] Outlier detection completed (Not removed).")


# =========================================================
#           Exploratory Data Analysis (EDA) 
# =========================================================

# --------------------------
# PCA (log-expression dataset)
# --------------------------

pca_res <- prcomp(t(expr_log), center = TRUE, scale. = TRUE)

pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = group,
  subtype = subtype,
  outlier = colnames(expr_log) %in% outlier_samples
)

pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

# PCA coloured by disease group

pca_plot_group <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = outlier)) +
  geom_point(size = 2.2, alpha = 0.85) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  scale_color_manual(values = c("#CE2915", "#0096FF")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  ggtitle("PCA of GSE19429 samples") +
  theme_bw()

# PCA colored by FAB subtype

n_subtypes <- nlevels(pca_df$subtype)
subtype_palette <- setNames(
  colorRampPalette(brewer.pal(min(max(n_subtypes, 3), 9), "Set2"))(n_subtypes),
  levels(pca_df$subtype)
)
pca_plot_subtype <- ggplot(pca_df, aes(PC1, PC2, color = subtype, shape = outlier)) +
  geom_point(size = 2.2, alpha = 0.85) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  scale_color_manual(values = subtype_palette) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  ggtitle("PCA of GSE19429 samples, colored by FAB subtype") +
  theme_bw()

save_plot(pca_plot_group, filename = "03_pca_plot.png", dir = paths$plots_gse19429_eda)
save_plot(pca_plot_subtype, filename = "03.1_pca_plot_subtype.png", dir = paths$plots_gse19429_eda)

# --------------------------
# Scree plot
# --------------------------

var_df <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  variance = pca_var[1:10] * 100,
  cumvar = cumsum(pca_var[1:10]) * 100
)

open_png(filename = "04_pca_scree_plot.png", dir = paths$plots_gse19429_eda,
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

sample_dist <- dist(t(expr_log), method = "euclidean")

hc <- hclust(sample_dist, method = "average")

# plot

open_png(filename = "05_hierarchical_clustering.png", dir = paths$plots_gse19429_eda,
  width = 2200, height = 900)

par(mar = c(5, 4, 4, 2))

plot(hc, labels = FALSE,
  main = "Hierarchical clustering of samples (Euclidean, average linkage)",
  xlab = "",
  sub = ""
)

# Color annotation aligned to dendrogram order

ordered_groups <- group[hc$order]

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

sample_cor <- cor(expr_log, method = "pearson", use = "pairwise.complete.obs")

# Sample annotations

annotation_col <- data.frame(Group = group, Subtype = subtype)
rownames(annotation_col) <- colnames(expr_log)

# Ensure annotation order matches correlation matrix

annotation_col <- annotation_col[colnames(sample_cor), , drop = FALSE]

stopifnot(identical(rownames(annotation_col), colnames(cor_mat)))
stopifnot(all(rownames(annotation_col) == colnames(cor_mat)))

open_png(filename = "06_sample_correlation_heatmap.png", dir = paths$plots_gse19429_eda,
         width = 1400, height = 1400)

pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_palette, Subtype = subtype_palette),
  show_rownames = FALSE, 
  show_colnames = FALSE,
  main = "Sample-sample Pearson correlation (GSE19429)"
)

close_png()


# ---------------------------------
# Sample-distance heatmap (1 - correlation)
# ---------------------------------

message("[QC] Computing sample distance matrix (1 - cor)...")

sample_dist_1mcor <- as.dist(1 - sample_cor)

open_png(filename = "07_sample_distance_heatmap.png", dir = paths$plots_gse19429_eda,
         width = 1400, height = 1400)
pheatmap(
  as.matrix(sample_dist_1mcor),
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_palette, Subtype = subtype_palette),
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "Sample distance heatmap (1 - Pearson correlation)"
)

close_png()


# =========================================================
#           Low-variance probe filtering 
# =========================================================

message("[QC] Starting low-variance probe filtering...")

# ---------------------------------------------------------
# Compute variability metric (IQR per probe)
# ---------------------------------------------------------

message("[QC] Computing probe-wise interquartile ranges...")

probe_iqr <- matrixStats::rowIQRs(expr_log, na.rm = TRUE)

iqr_cutoff <- quantile(probe_iqr, probs = 0.25, na.rm = TRUE)

keep_probes <- probe_iqr > iqr_cutoff

# ---------------------------------------------------------
# Reporting
# ---------------------------------------------------------

message(sprintf(
  "[QC] Low-variance filtering applied: keeping %d / %d probes (IQR cutoff = %.4f).",
  sum(keep_probes),
  length(keep_probes),
  iqr_cutoff
))

message(sprintf(
  "[QC] Removed %d low-variance probes.",
  sum(!keep_probes)
))

# ---------------------------------------------------------
# Generate filtered expression matrix
# ---------------------------------------------------------

expr_filt <- expr_log[keep_probes, , drop = FALSE]

message("[QC] Expression matrix filtered successfully.")


# =========================================================
#           DIFFERENTIAL EXPRESSION ANALYSIS (DEGs)
# =========================================================

stopifnot(all(colnames(expr_filt) == rownames(meta)))

message("[DEG] Building design matrix...")

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

stopifnot(all(c("MDS", "HC") %in% colnames(design)))

message("[DEG] Fitting linear model with limma...")

fit <- lmFit(expr_filt, design)

message("[DEG] Building contrast matrix...")

contrast_matrix <- makeContrasts(
  MDS_vs_HC = MDS - HC, 
  levels = design
  )

fit2 <- contrasts.fit(fit, contrast_matrix)

message("[DEG] Empirical Bayes moderation...")

fit2 <- eBayes(fit2, robust = TRUE)  

message("[DEG] Extracting DEG table...")

deg_table <- topTable(
  fit2, 
  coef = "MDS_vs_HC", 
  number = Inf, 
  adjust.method = "BH"
  )

deg_table$probe_id <- rownames(deg_table)

message(sprintf(
  "[DEG] Completed: %d features tested",
  nrow(deg_table)
))


# =========================================================
#            Probe annotation and DEG summarization
# =========================================================

message("[DEG] Annotating probes with gene information...")

entrez_symbol <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = deg_table$probe_id,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "PROBEID"
)

message("[DEG] Removing probes without gene annotation...")

entrez_symbol <- entrez_symbol[
  !(is.na(entrez_symbol$SYMBOL) &
      is.na(entrez_symbol$ENTREZID) &
      is.na(entrez_symbol$GENENAME)),
]

# ------------------------------------------
# Remove duplicated probe annotations (one annotation per probe)
# ------------------------------------------

message("[DEG] Removing probes with ambiguous gene annotations...")

entrez_symbol <- entrez_symbol[!duplicated(entrez_symbol$PROBEID), ]

message(sprintf(
  "[DEG] Retrieved annotations for %d unique probes.",
  nrow(entrez_symbol)
))

# ------------------------------------------
# Merge annotations with DEG results
# ------------------------------------------

deg_table <- merge(
  deg_table,
  entrez_symbol,
  by.x = "probe_id",
  by.y = "PROBEID",
  all.x = TRUE
)

deg_table <- deg_table[order(deg_table$adj.P.Val), ]

message("[DEG] Probe annotation completed.")

# ------------------------------------------
# Significant DEG summary
# ------------------------------------------

logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg <- subset(
  deg_table,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

n_before <- nrow(sig_deg)

sig_deg <- sig_deg[!(is.na(sig_deg$SYMBOL) & is.na(sig_deg$ENTREZID) & is.na(sig_deg$GENENAME)), ]

message(sprintf(
  "[DEG] Removed %d significant probe(s) with no annotation at all (SYMBOL/ENTREZID/GENENAME all NA).",
  n_before - nrow(sig_deg)
))

message(sprintf(
  "[DEG] Significant genes (FDR < %.2f, |log2FC| > %.1f): %d upregulated, %d downregulated (%d total).",
  padj_cutoff,
  logcf_cutoff,
  sum(sig_deg$logFC > 0),
  sum(sig_deg$logFC < 0),
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
deg_table$significance[deg_table$adj.P.Val < padj_cutoff & deg_table$logFC > logcf_cutoff] <- "Up in MDS"
deg_table$significance[deg_table$adj.P.Val < padj_cutoff & deg_table$logFC < -logcf_cutoff] <- "Down in MDS"

volcano_plot_deg <- ggplot(deg_table,aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
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
    title = "Volcano plot of differential expression (GSE19429)",
    x = expression(log[2] ~ fold~change),
    y = expression(-log[10] ~ adjusted~italic(P)),
    color = NULL
  ) +
  theme_bw()

save_plot(volcano_plot_deg, filename = "08_volcano_plot.png", dir = paths$plots_gse19429_deg, 
          width = 8, height = 6)

# ------------------------------------------
# MA Plot
# ------------------------------------------

message("[DEG] Generating MA plot...")

open_png(filename = "09_MA_plot_MDS_vs_HC.png", dir = paths$plots_gse19429_deg,
         width = 1000, height = 800)

limma::plotMA(fit2, coef = "MDS_vs_HC", status = deg_table$significance[match(rownames(fit2), deg_table$probe_id)],
              values = c("Up in MDS", "Down in MDS"), col = c("firebrick", "steelblue"),
              main = "MA plot: MDS vs Healthy Controls")
abline(h = 0, col = "grey40", lty = 2)

close_png()

# -------------------
# Heatmap
# ------------------

# Select top DEGs

top_degs <- 50

# Filter out non-annotated probes and re-order by statistical significance

deg_annotated <- sig_deg[!is.na(sig_deg$SYMBOL) & sig_deg$SYMBOL != "", ]
deg_annotated <- deg_annotated[order(deg_annotated$adj.P.Val), ]

# Deduplicate by gene symbol: multiple probes can map to the same SYMBOL

deg_annotated <- deg_annotated[!duplicated(deg_annotated$SYMBOL), ]

message(sprintf(
  "[DEG] %d unique annotated genes available after probe-to-gene deduplication.",
  nrow(deg_annotated)
))

# Sanity check

if (nrow(deg_annotated) < top_degs) {
  warning(sprintf(
    "[DEG] Only %d significant annotated genes available (requested %d) - heatmap will show all %d.",
    nrow(deg_annotated), top_degs, nrow(deg_annotated)
  ))
  top_degs <- nrow(deg_annotated)
}

message(sprintf(
  "[DEG] %d unique significant genes selected for heatmap (out of %d significant DEGs total).",
  top_degs, nrow(sig_deg)
))

# Final selection for heatmap:
# - uses only annotated probes
# - retains top statistically significant features

top_probes <- head(deg_annotated$probe_id, top_degs)
top_symbols <- head(deg_annotated$SYMBOL, top_degs)

heatmap_matrix <- expr_filt[top_probes, , drop = FALSE]

# Row-wise z-score scaling
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))

# Ensure unique gene symbols
rownames(heatmap_matrix_scaled) <- make.unique(top_symbols)

# Sample annotations
ann_col <- HeatmapAnnotation(
  Group = group,
  col = list(Group = group_palette)
)

# Plot
open_png(filename = "10-Heatmap_TopDegs_MDS_vs_HC.png", dir = paths$plots_gse19429_deg,
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

message("[OUTPUT] Saving differential expression results...")

save_csv(deg_table, filename = "DEGs_GSE19429.csv", dir = paths$tables)

message("[OUTPUT] Saving background gene list for enrichment analyses...")

background <- rownames(expr_filt)

save_csv(background, filename = "background_GSE19429.csv", dir = paths$tables)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_DEG_GSE19429.txt", dir = paths$logs, label = "Differential expression analysis - GSE19429")

message("[OUTPUT] Session information saved to:" , paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Differential expression analysis completed successfully for GSE19429 (MDS vs Healthy Controls).")
message("=================================================")
