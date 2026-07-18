# =========================================================
# Project      : MDS Genomics - Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE114922 (Illumina HiSeq 4000)
# Samples      : 82 MDS patients + 8 healthy controls, CD34+ bone marrow HSC
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

message("=== Starting GSE114922 RNA-seq preprocessing and QC ===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Loading raw count matrix...")

counts_raw <- read.delim(gzfile("Raw_data/GSE114922_raw_counts_GRCh38.p13_NCBI.tsv.gz"), row.names = 1)

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
#     SVA DIAGNOSTIC CORRECTION (RNA-seq) - GSE114922
# Methods: Buja and Eyuboglu Method and Leek
# =========================================================

message("[SVA] Extracting normalized counts...")

dat <- counts(dds, normalized = TRUE)
keep_nonzero <- rowSums(dat) > 0
dat <- dat[keep_nonzero, ]

message(sprintf(
  "[SVA] %d / %d features retained for SVA (non-zero rows).",
  sum(keep_nonzero), length(keep_nonzero)
))

mod  <- model.matrix(~ condition, data = colData(dds))
mod0 <- model.matrix(~ 1, data = colData(dds))

stopifnot(nrow(mod) == ncol(dat))

message("[SVA] Estimating number of surrogate variables (BE vs Leek)...")

n_sv_be   <- num.sv(dat, mod, method = "be")
n_sv_leek <- num.sv(dat, mod, method = "leek")

message(sprintf("[SVA] n.sv (be): %d | n.sv (leek): %d", n_sv_be, n_sv_leek))

set.seed(1234)

svseq_be <- tryCatch({
  if (n_sv_be > 0) {
    svaseq(dat, mod, mod0, n.sv = n_sv_be)
  } else {
    NULL
  }
}, error = function(e) {
  message(sprintf(
    "[SVA] Leek method failed to converge (n.sv=%d requested on %d samples): %s",
    n_sv_leek, ncol(dat), conditionMessage(e)
  ))
  NULL
})

svseq_leek <- tryCatch({
  if (n_sv_leek > 0) {
    svaseq(dat, mod, mod0, n.sv = n_sv_leek)
  } else {
    NULL
  }
}, error = function(e) {
  message(sprintf(
    "[SVA] Leek method failed to converge (n.sv=%d requested on %d samples): %s",
    n_sv_leek, ncol(dat), conditionMessage(e)
  ))
  NULL
})


# -----------
# PCA plot
# -----------

make_sva_pca <- function(mat, covariates, title_suffix) {
  
  mat_corrected <- removeBatchEffect(mat, covariates = covariates, design = mod)
  
  pca_res <- prcomp(t(mat_corrected), center = TRUE, scale. = FALSE)
  pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
  pca_var_percent <- round(100 * pca_var / sum(pca_var), 1)
  
  pca_df_sva <- data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    group = meta_cd34$group
  )
  
  ggplot(pca_df_sva, aes(PC1, PC2, color = group)) +
    geom_point(size = 2.2, alpha = 0.85) +
    xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
    ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
    scale_color_manual(values = c("#CE2915", "#0096FF")) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
    ggtitle(paste0("PCA of GSE114922 samples - SVA (", title_suffix, ")")) +
    theme_bw()
}

if (!is.null(svseq_be)) {
  pca_plot_sva_be <- make_sva_pca(vsdata_df, svseq_be$sv, sprintf("Buja and Eyuboglu Method, SV = %d", svseq_be$n.sv))
  save_plot(pca_plot_sva_be, filename = "pca_plot_SVA_be.png", dir = paths$gse_114922_sva)
}

if (!is.null(svseq_leek)) {
  pca_plot_sva_leek <- make_sva_pca(vsdata_df, svseq_leek$sv, sprintf("Leek, n.sv=%d", svseq_leek$n.sv))
  save_plot(pca_plot_sva_leek, filename = "pca_plot_SVA_leek.png", dir = paths$plots_gse114922_eda)
} else message("NULL")


# -----------
# DEG - BE 
# -----------

logcf_cutoff <- 1
padj_cutoff  <- 0.05

deg_table_sva_be <- NULL
sig_deg_sva_be   <- NULL

if (!is.null(svseq_be)) {
  
  message("[DEG-SVA] Fitting model with BE surrogate variables...")
  
  dds_be <- dds
  for (i in seq_len(svseq_be$n.sv)) {
    colData(dds_be)[[paste0("SV", i)]] <- svseq_be$sv[, i]
  }
  sv_terms_be <- paste0("SV", seq_len(svseq_be$n.sv), collapse = " + ")
  design(dds_be) <- as.formula(paste("~", sv_terms_be, "+ condition"))
  
  dds_be <- DESeq(dds_be)
  
  deg_table_sva_be_raw <- results(dds_be, contrast = c("condition", "MDS", "HC"), alpha = 0.05)
  deg_table_sva_be <- as.data.frame(deg_table_sva_be_raw)
  deg_table_sva_be$ENTREZID <- rownames(deg_table_sva_be)
  
  deg_table_sva_be$SYMBOL <- mapIds(
    org.Hs.eg.db, keys = rownames(deg_table_sva_be),
    keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"
  )
  deg_table_sva_be$GENENAME <- mapIds(
    org.Hs.eg.db, keys = rownames(deg_table_sva_be),
    keytype = "ENTREZID", column = "GENENAME", multiVals = "first"
  )
  
  sig_deg_sva_be <- subset(
    deg_table_sva_be,
    padj < padj_cutoff & abs(log2FoldChange) > logcf_cutoff
  )
  sig_deg_sva_be <- sig_deg_sva_be[!is.na(sig_deg_sva_be$padj), ]
  sig_deg_sva_be <- sig_deg_sva_be[
    !(is.na(sig_deg_sva_be$SYMBOL) & is.na(sig_deg_sva_be$ENTREZID) & is.na(sig_deg_sva_be$GENENAME)), 
  ]
  
  message(sprintf(
    "[DEG-SVA] BE: %d upregulated, %d downregulated (%d total).",
    sum(sig_deg_sva_be$log2FoldChange > 0),
    sum(sig_deg_sva_be$log2FoldChange < 0),
    nrow(sig_deg_sva_be)
  ))
}


# -----------
# DEG - Leek
# -----------

deg_table_sva_leek <- NULL
sig_deg_sva_leek   <- NULL

if (!is.null(svseq_leek)) {
  
  message("[DEG-SVA] Fitting model with Leek surrogate variables...")
  
  dds_leek <- dds
  for (i in seq_len(svseq_leek$n.sv)) {
    colData(dds_leek)[[paste0("SV", i)]] <- svseq_leek$sv[, i]
  }
  sv_terms_leek <- paste0("SV", seq_len(svseq_leek$n.sv), collapse = " + ")
  design(dds_leek) <- as.formula(paste("~", sv_terms_leek, "+ condition"))
  
  dds_leek <- DESeq(dds_leek)
  
  deg_table_sva_leek_raw <- results(dds_leek, contrast = c("condition", "MDS", "HC"), alpha = 0.05)
  deg_table_sva_leek <- as.data.frame(deg_table_sva_leek_raw)
  deg_table_sva_leek$ENTREZID <- rownames(deg_table_sva_leek)
  
  deg_table_sva_leek$SYMBOL <- mapIds(
    org.Hs.eg.db, keys = rownames(deg_table_sva_leek),
    keytype = "ENTREZID", column = "SYMBOL", multiVals = "first"
  )
  deg_table_sva_leek$GENENAME <- mapIds(
    org.Hs.eg.db, keys = rownames(deg_table_sva_leek),
    keytype = "ENTREZID", column = "GENENAME", multiVals = "first"
  )
  
  sig_deg_sva_leek <- subset(
    deg_table_sva_leek,
    padj < padj_cutoff & abs(log2FoldChange) > logcf_cutoff
  )
  sig_deg_sva_leek <- sig_deg_sva_leek[!is.na(sig_deg_sva_leek$padj), ]
  sig_deg_sva_leek <- sig_deg_sva_leek[
    !(is.na(sig_deg_sva_leek$SYMBOL) & is.na(sig_deg_sva_leek$ENTREZID) & is.na(sig_deg_sva_leek$GENENAME)), 
  ]
  
  message(sprintf(
    "[DEG-SVA] Leek: %d upregulated, %d downregulated (%d total).",
    sum(sig_deg_sva_leek$log2FoldChange > 0),
    sum(sig_deg_sva_leek$log2FoldChange < 0),
    nrow(sig_deg_sva_leek)
  ))
} else message("NULL")


# -------------
# Summary
# -------------

n_deg_none <- sum(res$padj < 0.05 & abs(res$log2FoldChange) > 1, na.rm = TRUE)

message(sprintf(
  "[SVA] DEG comparison - None: %d | BE (n.sv=%d): %s | Leek (n.sv=%d): %s",
  n_deg_none,
  n_sv_be,   ifelse(is.null(sig_deg_sva_be),   "skipped", nrow(sig_deg_sva_be)),
  n_sv_leek, ifelse(is.null(sig_deg_sva_leek), "skipped", nrow(sig_deg_sva_leek))
))


# =========================================================
#                Save outputs
# =========================================================

message("[OUTPUT] Saving differential expression results...")

save_csv(deg_table_sva_be, filename = "DEGs_SVA_BE_GSE114922.csv", dir = paths$tables_sva)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_DEG_SVA_GSE114922.txt", dir = paths$logs_sva, label = "Build DESeq2 object - SVA Correction - GSE114922")

message("[OUTPUT] Session information saved to:" , paths$logs_sva)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Differential expression analysis, using SVA correction, completed successfully for GSE114922 (MDS vs Healthy Controls).")
message("=================================================")

