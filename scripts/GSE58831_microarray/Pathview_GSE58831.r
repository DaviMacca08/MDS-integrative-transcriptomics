# =========================================================
# Project      : MDS Genomics - Microarray Analysis Pipeline
# Dataset      : GSE58831 (Affymetrix HG-U133 Plus 2.0, GPL570)
# Samples      : 183 MDS patients + 17 healthy controls, CD34+ bone marrow HSC
# Script       : Pathway Visualization with Pathview
# Description  : Visualization of selected biological pathways by overlaying
#                gene-level expression changes (log2 fold changes) onto KEGG
#                pathway maps using Pathview.
#                Pathway selection is based on previous Differential Expression
#                Analysis (DEG), Over-Representation Analysis (ORA), and
#                Gene Set Enrichment Analysis (GSEA) results.
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("SetupEnvironment/00_paths.R")
source("SetupEnvironment/01_environment.R")
source("SetupEnvironment/02_seed.R")
source("SetupEnvironment/03_helper_functions.R")

set_seed(1234)

message("=== Starting GSE58831 Pathview===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Reading differential expression results...")

deg_table <- read.csv(file = file.path(paths$tables, "DEGs_GSE58831.csv"), sep = "," )

# Prepare data for pathview

list_pathview <- deg_table$logFC
names(list_pathview) <- deg_table$ENTREZID


# =========================================================
#                  Visualize pathway 
# =========================================================
# JAK-STAT signaling

withr::with_dir(
  "results/GSE58831_microarray/PathView/",
  {
    pathview(
      gene.data = list_pathview,
      pathway.id = "hsa04630",
      species = "hsa",
      res = 600,
      cex = 0.4,
      low  = list(gene = "blue"),
      mid  = list(gene = "white"),
      high = list(gene = "red"),
      new.signature = FALSE,
      out.suffix = "JAK_STAT_MDS_vs_HC"
    )
  }
)

# NF-kB signaling

withr::with_dir(
  "results/GSE58831_microarray/PathView/",
  {
    pathview(
      gene.data = list_pathview,
      pathway.id = "hsa04064",
      species = "hsa",
      res = 600,
      cex = 0.4,
      low  = list(gene = "blue"),
      mid  = list(gene = "white"),
      high = list(gene = "red"),
      new.signature = FALSE,
      out.suffix = "NFkB_MDS_vs_HC"
    )
  }
)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_Pathview_GSE58831.txt", dir = paths$logs, label = "Pathway Visualization with Pathview - GSE58831")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Pathway visualization with Pathview completed successfully for GSE58831 (MDS vs Healthy Controls).")
message("=================================================")

