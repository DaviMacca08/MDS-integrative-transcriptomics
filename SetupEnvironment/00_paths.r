# =========================================================
# Project paths (single source of truth)
# =========================================================

base_dir <- "/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects🤞🏻/Expression-profiles-MDS/github"

paths <- list(
  
  base = base_dir,

  # results split by dataset
  results = file.path(base_dir, "results"),
  gse_19429 = file.path(base_dir, "results/GSE19429"),
  gse_114922 = file.path(base_dir, "results/GSE114922"),
  gse_58831 = file.path(base_dir, "results/GSE58831"),
  
  # plots split by dataset
  plots_gse19429_qc = file.path(base_dir, "results/GSE19429/QC"),
  plots_gse19429_eda = file.path(base_dir, "results/GSE19429/EDA"),
  plots_gse19429_deg = file.path(base_dir, "results/GSE19429/DEG"),
  plots_gse19429_ora = file.path(base_dir, "results/GSE19429/ORA"),
  plots_gse19429_gsea = file.path(base_dir, "results/GSE19429/GSEA"),
  plots_gse19429_path = file.path(base_dir, "results/GSE19429/PathView"),
  
  plots_gse114922_qc = file.path(base_dir, "results/GSE114922/QC"),
  plots_gse114922_eda = file.path(base_dir, "results/GSE114922/EDA"),
  plots_gse114922_deg = file.path(base_dir, "results/GSE114922/DEG"),
  plots_gse114922_ora = file.path(base_dir, "results/GSE114922/ORA"),
  plots_gse114922_gsea = file.path(base_dir, "results/GSE114922/GSEA"),
  plots_gse114922_path = file.path(base_dir, "results/GSE114922/PathView"),
  
  plots_gse58831_qc = file.path(base_dir, "results/GSE58831/QC"),
  plots_gse58831_eda = file.path(base_dir, "results/GSE58831/EDA"),
  plots_gse58831_deg = file.path(base_dir, "results/GSE58831/DEG"),
  plots_gse58831_ora = file.path(base_dir, "results/GSE58831/ORA"),
  plots_gse58831_gsea = file.path(base_dir, "results/GSE58831/GSEA"),
  plots_gse58831_path = file.path(base_dir, "results/GSE58831/PathView"),
  
  # tables inside results
  tables = file.path(base_dir, "results/tables"),
  
  # logs inside results
  logs = file.path(base_dir, "results/logs")
  
)

# Ensure directories exist
invisible(lapply(paths, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}))
