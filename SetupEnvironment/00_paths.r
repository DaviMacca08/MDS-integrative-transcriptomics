# =========================================================
# Project paths (single source of truth)
# =========================================================

base_dir <- "/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects🤞🏻/MDS-integrative-transcriptomics"

paths <- list(
  
  base = base_dir,

#------------------
# No SVA correction
#------------------

  # results split by dataset
  results = file.path(base_dir, "results"),
  gse_19429 = file.path(base_dir, "results/GSE19429_microarray"),
  gse_114922 = file.path(base_dir, "results/GSE114922_rnaseq"),
  gse_58831 = file.path(base_dir, "results/GSE58831_microarray"),
  
  # plots split by dataset
  plots_gse19429_qc = file.path(base_dir, "results/GSE19429_microarray/QC"),
  plots_gse19429_eda = file.path(base_dir, "results/GSE19429_microarray/EDA"),
  plots_gse19429_deg = file.path(base_dir, "results/GSE19429_microarray/DEG"),
  plots_gse19429_ora = file.path(base_dir, "results/GSE19429_microarray/ORA"),
  plots_gse19429_gsea = file.path(base_dir, "results/GSE19429_microarray/GSEA"),
  plots_gse19429_path = file.path(base_dir, "results/GSE19429_microarray/PathView"),
  
  plots_gse114922_qc = file.path(base_dir, "results/GSE114922_rnaseq/QC"),
  plots_gse114922_eda = file.path(base_dir, "results/GSE114922_rnaseq/EDA"),
  plots_gse114922_deg = file.path(base_dir, "results/GSE114922_rnaseq/DEG"),
  plots_gse114922_ora = file.path(base_dir, "results/GSE114922_rnaseq/ORA"),
  plots_gse114922_gsea = file.path(base_dir, "results/GSE114922_rnaseq/GSEA"),
  plots_gse114922_path = file.path(base_dir, "results/GSE114922_rnaseq/PathView"),
  
  plots_gse58831_qc = file.path(base_dir, "results/GSE58831_microarray/QC"),
  plots_gse58831_eda = file.path(base_dir, "results/GSE58831_microarray/EDA"),
  plots_gse58831_deg = file.path(base_dir, "results/GSE58831_microarray/DEG"),
  plots_gse58831_ora = file.path(base_dir, "results/GSE58831_microarray/ORA"),
  plots_gse58831_gsea = file.path(base_dir, "results/GSE58831_microarray/GSEA"),
  plots_gse58831_path = file.path(base_dir, "results/GSE58831_microarray/PathView"),
  
  # plots comparison
  plots_comparison = file.path(base_dir, "results/Comparison/"),
  
  # tables inside results
  tables = file.path(base_dir, "results/tables"),
  tables_gsea = file.path(base_dir, "results/tables/GSEA_results"),
  
  # logs inside results
  logs = file.path(base_dir, "results/logs"),
  
#------------------
# SVA correction
#------------------

# results split by dataset
results_sva = file.path(base_dir, "results_sva"),
gse_19429_sva = file.path(base_dir, "results_sva/GSE19429_microarray"),
gse_114922_sva = file.path(base_dir, "results_sva/GSE114922_rnaseq"),
gse_58831_sva = file.path(base_dir, "results_sva/GSE58831_microarray"),

# tables inside results
tables_sva = file.path(base_dir, "results_sva/tables"),

# logs inside results
logs_sva = file.path(base_dir, "results_sva/logs")


)

# Ensure directories exist
invisible(lapply(paths, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}))
