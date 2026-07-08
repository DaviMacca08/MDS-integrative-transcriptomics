# =========================================================
# Environment & packages
# =========================================================

suppressPackageStartupMessages({
  
  # Data manipulation & visualization
  library(tidyverse)
  library(matrixStats)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
  
  # GEO data acquisition & annotation
  library(GEOquery)
  library(hgu133plus2.db)
  
  # Microarray differential expression
  library(limma)
  
  # RNA-seq differential expression
  library(DESeq2)
  
  # Heatmaps & visualization
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  
  # Functional enrichment analysis
  library(clusterProfiler)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(fgsea)
  
  # Pathway visualization
  library(pathview)
  
})