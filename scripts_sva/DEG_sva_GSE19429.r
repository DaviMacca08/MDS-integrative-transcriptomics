# =========================================================
# Project      : MDS Genomics - Microarray Analysis Pipeline
# Dataset      : GSE19429 (Affymetrix HG-U133 Plus 2.0, GPL570)
# Samples      : 183 MDS patients + 17 healthy controls, CD34+ bone marrow HSC
# Script       : Differential Expression Analysis after SVA Correction
# Description  : Differential expression analysis (MDS vs Healthy Controls)
#                following surrogate variable analysis (SVA) using the
#                Buja & Eyuboglu (BE) and Leek methods.
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
#           SVA Correction
# =========================================================

# Ensure expression matrix and metadata are properly aligned
stopifnot(identical(colnames(expr_filt), rownames(meta)))

# Full model including the biological variable of interest
mod <- model.matrix(~ group, data = meta)
stopifnot(nrow(mod) == ncol(expr_filt))

# Null model (intercept only)
mod0 <- model.matrix(~ 1, data = meta)

# Estimate the number of surrogate variables using the Leek method
n_sv_leek <- num.sv(expr_filt, mod, method = "leek") 
n_sv_be <- num.sv(expr_filt, mod, method = "be")

message("[SVA] Estimated number of surrogate variables (Be method): ", n_sv_be)
message("[SVA] Estimated number of surrogate variables (Leek method): ", n_sv_leek)

# Estimate surrogate variables
sv_leek <- sva(expr_filt, mod, mod0, n.sv = n_sv_leek)
sv_be <- sva(expr_filt, mod, mod0, n.sv = n_sv_be)

message("[SVA] SVA completed successfully.")
message("[SVA] Number of surrogate variables estimated: ", sv_leek$n.sv)

# Remove surrogate variable effects for exploratory data visualization only
expr_sva_leek_vis <- removeBatchEffect(
  expr_filt,
  covariates = sv_leek$sv,
  design = mod
)

expr_sva_be_vis <- removeBatchEffect(
  expr_filt,
  covariates = sv_be$sv,
  design = mod
)


# =========================================================
#           Exploratory Data Analysis (EDA) 
# =========================================================

# ---------------------------------------------
# PCA (log-expression dataset) - SVA BE Method
# ---------------------------------------------

pca_res <- prcomp(t(expr_sva_be_vis), center = TRUE, scale. = TRUE)

pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = group
)

pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

# PCA coloured by disease group 

pca_plot_group_be <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  scale_color_manual(values = c("#CE2915", "#0096FF")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  ggtitle(paste0("PCA of GSE19429 Samples - SVA (Buja and Eyuboglu Method, SV = ", n_sv_be, ")")) +
  theme_bw()

save_plot(pca_plot_group_be, filename = "pca_plot_SVA_be.png", dir = paths$gse_19429_sva)

# ---------------------------------------------
# PCA (log-expression dataset) - SVA Leek Method
# ---------------------------------------------

pca_res <- prcomp(t(expr_sva_leek_vis), center = TRUE, scale. = TRUE)

pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = group
)

pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)

# PCA coloured by disease group 

pca_plot_group_leek <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2.2, alpha = 0.85) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  scale_color_manual(values = c("#CE2915", "#0096FF")) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  ggtitle(paste0("PCA of GSE19429 Samples - SVA (Leek Method SV = ", n_sv_leek, ")")) +
  theme_bw()

save_plot(pca_plot_group_leek, filename = "pca_plot_SVA_leek.png", dir = paths$gse_19429_sva)


# =========================================================
#           DIFFERENTIAL EXPRESSION ANALYSIS (DEGs)
#
# SVA Corrections using BE Method
# =========================================================

stopifnot(all(colnames(expr_filt) == rownames(meta)))

message("[DEG] Building design matrix...")

design <- model.matrix(~group + sv_be$sv)
colnames(design)[-(1:2)] <- paste0("SV", seq_len(sv_be$n.sv))

message("[DEG] Fitting linear model with limma...")

fit <- lmFit(expr_filt, design)

message("[DEG] Building contrast matrix...")

contrast_matrix <- makeContrasts(
  MDS_vs_HC = groupMDS, 
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

deg_table_sva_be <- merge(
  deg_table,
  entrez_symbol,
  by.x = "probe_id",
  by.y = "PROBEID",
  all.x = TRUE
)

deg_table_sva_be <- deg_table_sva_be[order(deg_table_sva_be$adj.P.Val), ]

message("[DEG] Probe annotation completed.")

# ------------------------------------------
# Significant DEG summary
# ------------------------------------------

logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg_sva_be <- subset(
  deg_table_sva_be,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

n_before <- nrow(sig_deg_sva_be)

sig_deg_sva_be <- sig_deg_sva_be[!(is.na(sig_deg_sva_be$SYMBOL) & is.na(sig_deg_sva_be$ENTREZID) & is.na(sig_deg_sva_be$GENENAME)), ]

message(sprintf(
  "[DEG] Removed %d significant probe(s) with no annotation at all (SYMBOL/ENTREZID/GENENAME all NA).",
  n_before - nrow(sig_deg_sva_be)
))

message(sprintf(
  "[DEG] Significant genes (FDR < %.2f, |log2FC| > %.1f): %d upregulated, %d downregulated (%d total).",
  padj_cutoff,
  logcf_cutoff,
  sum(sig_deg_sva_be$logFC > 0),
  sum(sig_deg_sva_be$logFC < 0),
  nrow(sig_deg_sva_be)
))


# =========================================================
#           DIFFERENTIAL EXPRESSION ANALYSIS (DEGs)
#
# SVA Corrections using LEEK Method
# =========================================================

stopifnot(all(colnames(expr_filt) == rownames(meta)))

message("[DEG] Building design matrix...")

design <- model.matrix(~group + sv_leek$sv)
colnames(design)[-(1:2)] <- paste0("SV", seq_len(sv_leek$n.sv))

message("[DEG] Fitting linear model with limma...")

fit <- lmFit(expr_filt, design)

message("[DEG] Building contrast matrix...")

contrast_matrix <- makeContrasts(
  MDS_vs_HC = groupMDS, 
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

deg_table_sva_leek <- merge(
  deg_table,
  entrez_symbol,
  by.x = "probe_id",
  by.y = "PROBEID",
  all.x = TRUE
)

deg_table_sva_leek <- deg_table_sva_leek[order(deg_table_sva_leek$adj.P.Val), ]

message("[DEG] Probe annotation completed.")

# ------------------------------------------
# Significant DEG summary
# ------------------------------------------

logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg_sva_leek <- subset(
  deg_table_sva_leek,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

n_before <- nrow(sig_deg_sva_leek)

sig_deg_sva_leek <- sig_deg_sva_leek[!(is.na(sig_deg_sva_leek$SYMBOL) & is.na(sig_deg_sva_leek$ENTREZID) & is.na(sig_deg_sva_leek$GENENAME)), ]

message(sprintf(
  "[DEG] Removed %d significant probe(s) with no annotation at all (SYMBOL/ENTREZID/GENENAME all NA).",
  n_before - nrow(sig_deg_sva_leek)
))

message(sprintf(
  "[DEG] Significant genes (FDR < %.2f, |log2FC| > %.1f): %d upregulated, %d downregulated (%d total).",
  padj_cutoff,
  logcf_cutoff,
  sum(sig_deg_sva_leek$logFC > 0),
  sum(sig_deg_sva_leek$logFC < 0),
  nrow(sig_deg_sva_leek)
))


# =========================================================
#                Save outputs
# =========================================================

message("[OUTPUT] Saving differential expression results...")

save_csv(deg_table_sva_be, filename = "DEGs_SVA_BE_GSE19429.csv", dir = paths$tables_sva)
save_csv(deg_table_sva_leek, filename = "DEGs_SVA_LEEK_GSE19429.csv", dir = paths$tables_sva)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_DEG_SVA_GSE19429.txt", dir = paths$logs_sva, label = "Differential expression analysis - SVA Correction - GSE19429")

message("[OUTPUT] Session information saved to:" , paths$logs_sva)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Differential expression analysis, using SVA correction, completed successfully for GSE19429 (MDS vs Healthy Controls).")
message("=================================================")

