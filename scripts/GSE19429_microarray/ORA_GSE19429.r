# =========================================================
# Project      : MDS Genomics - Microarray Analysis Pipeline
# Dataset      : GSE19429 (Affymetrix HG-U133 Plus 2.0, GPL570)
# Samples      : 183 MDS patients + 17 healthy controls, CD34+ bone marrow HSC
# Script       : Over-Representation Analysis (ORA)
# Description  : Functional over-representation analysis using GO, KEGG, Reactome and MSigDB Hallmark
# =========================================================

# =========================================================
#                     Functions
# =========================================================

# Create enrichment dot plot

plot_enrichment <- function(df, title, color, top_n = 10, wrap_width = 35, fdr_cutoff = 0.05) {
  
  if (is.null(df) || nrow(df) == 0) {
    message("[PLOT] Empty or NULL enrichment result. No plot generated.")
    return(NULL)
  }
  
  df <- as.data.frame(df)
  df <- df[!is.na(df$p.adjust), ]
  
  n_before <- nrow(df)
  df <- df[!is.na(df$p.adjust) & df$p.adjust < fdr_cutoff, ]
  
  message(sprintf(
    "[PLOT] %d / %d terms retained after FDR filter (p.adjust < %.2f).",
    nrow(df), n_before, fdr_cutoff
  ))
  
  if (nrow(df) == 0) {
    message("[PLOT] No enriched terms passed the FDR threshold.)")
    return(NULL)
  }
  
  if (!"FoldEnrichment" %in% colnames(df)) {
    
    message("[PLOT] Calculating Fold Enrichment...")
    
    gene_ratio_num <- as.numeric(sub("/.*", "", df$GeneRatio))
    gene_ratio_den <- as.numeric(sub(".*/", "", df$GeneRatio))
    
    bg_ratio_num <- as.numeric(sub("/.*", "", df$BgRatio))
    bg_ratio_den <- as.numeric(sub(".*/", "", df$BgRatio))
    
    df$FoldEnrichment <- (gene_ratio_num / gene_ratio_den) / (bg_ratio_num / bg_ratio_den)
  }
  
  has_set <- "Set" %in% colnames(df)
  has_ontology <- "ONTOLOGY" %in% colnames(df)
  
  if (has_set) {
    df$Set <- factor(df$Set, levels = c("Module", "DEG", "Intersection"))
  } else {
    df$GeneRatioNum <- as.numeric(sub("/.*", "", df$GeneRatio)) /
      as.numeric(sub(".*/", "", df$GeneRatio))
  }
  
  df$Description <- stringr::str_wrap(df$Description, width = wrap_width)
  
  group_vars <- intersect(c("Set", "ONTOLOGY"), colnames(df))
  
  if (length(group_vars) > 0) {
    Top <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::arrange(p.adjust, .by_group = TRUE) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  } else {
    Top <- df |>
      dplyr::arrange(p.adjust) |>
      dplyr::slice_head(n = top_n)
  }
  
  Top$Description <- factor(
    Top$Description,
    levels = Top$Description[order(Top$p.adjust)]
  )
  
  x_var <- if (has_set) "Set" else "GeneRatioNum"
  
  stopifnot(x_var %in% colnames(Top))
  
  p <- ggplot(Top, aes(x = .data[[x_var]], y = Description)) +
    geom_point(aes(size = FoldEnrichment, color = p.adjust)) +
    scale_color_gradient(low = color, high = "#2C7BB6", labels = scales::scientific) +
    scale_size_continuous(range = c(2, 6), breaks = pretty(Top$FoldEnrichment, n = 3)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
      strip.text.y = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = title,
      x = if (has_set) "" else "Gene Ratio",
      y = "",
      size = "Fold Enrichment",
      color = "FDR"
    )
  
  if (has_ontology) {
    p <- p + facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y")
  }
  
  message("[PLOT] Top terms plotted: ", nrow(Top))
  
  p
  
}


# =========================================================
#                  Libraries & Setup
# =========================================================

source("SetupEnvironment/00_paths.R")
source("SetupEnvironment/01_environment.R")
source("SetupEnvironment/02_seed.R")
source("SetupEnvironment/03_helper_functions.R")

set_seed(1234)

message("=== Starting over-representation analysis for GSE19429 ===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Reading differential expression results...")

deg_table <- read.csv(file = file.path(paths$tables, "DEGs_GSE19429.csv"), sep = "," )

message("[LOAD] Reading background gene list...")

background <- read.csv(file = file.path(paths$tables, "background_GSE19429.csv"), sep = "," )


# =========================================================
#  Prepare DEG and background gene lists for enrichment analysis
# =========================================================

logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg <- subset(
  deg_table,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

n_before <- nrow(sig_deg)

sig_deg <- sig_deg[!(is.na(sig_deg$SYMBOL) & is.na(sig_deg$ENTREZID) & is.na(sig_deg$GENENAME)), ]

# Up-regultaed genes

sig_deg_up <- sig_deg |> 
  filter(logFC > 0)

# Down-regultaed genes

sig_deg_down <- sig_deg |> 
  filter(logFC < 0)

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

# Create Entrez ID gene lists

list_deg <- unique(sig_deg$ENTREZID[!is.na(sig_deg$ENTREZID)])
list_deg_up <- unique(sig_deg_up$ENTREZID[!is.na(sig_deg_up$ENTREZID)])
list_deg_down <- unique(sig_deg_down$ENTREZID[!is.na(sig_deg_down$ENTREZID)])

list_universe_raw <- AnnotationDbi::select(
  hgu133plus2.db,
  keys = background$x,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = "PROBEID"
)

list_universe <- unique(list_universe_raw$ENTREZID[!is.na(list_universe_raw$ENTREZID)])


# =========================================================
#                  Run ORA - GO Database
# =========================================================

# Gene Ontology over-representation analysis (BP, MF and CC)
# GO Enrichment - all DEGs

ego <- enrichGO(
  gene = list_deg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
  )

# GO Enrichment - upregulated DEGs

ego_up <- enrichGO(
  gene = list_deg_up,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# GO Enrichment - downregulated DEGs

ego_down <- enrichGO(
  gene = list_deg_down,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "ALL",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# Generate enrichment plots

go_all_plot <- plot_enrichment(ego, "GO Enrichment - MDS vs HC (All genes)", color = "#C44E52")
go_up_plot <- plot_enrichment(ego_up, "GO Enrichment - MDS vs HC (Up-regulated genes)", color = "#009E73")
go_down_plot <- plot_enrichment(ego_down, "GO Enrichment - MDS vs HC (Down-regulated genes)", color = "#CC79A7")

# Save enrichment plots

save_plot(go_all_plot, filename = "11-GO_all_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(go_up_plot, filename = "12-GO_Up_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(go_down_plot, filename = "13-GO_Down_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)


# =========================================================
#                  Run ORA - KEGG Database
# =========================================================

# KEGG enrichment - all DEGs

kegg_all <- enrichKEGG(
  gene = list_deg,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
  )

# KEGG enrichment - upregulated DEGs

kegg_up <- enrichKEGG(
  gene = list_deg_up,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# KEGG enrichment - downregulated DEGs

kegg_down <- enrichKEGG(
  gene = list_deg_down,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# Generate enrichment plots

kegg_all_plot <- plot_enrichment(kegg_all, title = "KEGG Enrichment - MDS vs HC (All genes)", color = "#C44E52")
kegg_up_plot <- plot_enrichment(kegg_up, title = "KEGG Enrichment - MDS vs HC (Up-regulated genes)", color = "#009E73")
kegg_down_plot <- plot_enrichment(kegg_down, title = "KEGG Enrichment - MDS vs HC (Down-regulated genes)", color = "#CC79A7")

# Save enrichment plots

save_plot(kegg_all_plot, filename = "14-KEGG_all_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(kegg_up_plot, filename = "15-KEGG_up_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(kegg_down_plot, filename = "16-KEGG_down_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)


# =========================================================
#                  Run ORA - Reactome Database
# =========================================================

# Reactome enrichment - all DEGs

react_all <- enrichPathway(
  gene = list_deg,
  organism = "human",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

# Reactome enrichment - upregulated DEGs

react_up <- enrichPathway(
  gene = list_deg_up,
  organism = "human",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

# Reactome enrichment - downregulated DEGs

react_down <- enrichPathway(
  gene = list_deg_down,
  organism = "human",
  universe = list_universe,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

# Generate enrichment plots

react_all_plot <- plot_enrichment(react_all, title = "Reactome Enrichment - MDS vs HC (All genes)", color = "#C44E52")
react_up_plot <- plot_enrichment(react_up, title = "Reactome Enrichment - MDS vs HC (Up-regulated genes)", color = "#009E73")
react_down_plot <- plot_enrichment(react_down, title = "Reactome Enrichment - MDS vs HC (Down-regulated genes)", color = "#CC79A7")

# Save enrichment plots

save_plot(react_all_plot, filename = "17-Reactome_all_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(react_up_plot, filename = "18-Reactome_up_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(react_down_plot, filename = "19-Reactome_down_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)


# =========================================================
#                  Run ORA - MSigDB Hallmark
# =========================================================

message("[ORA] Retrieving MSigDB Hallmark (H) gene sets...")

if ("collection" %in% names(formals(msigdbr::msigdbr))) {
  hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  message("Collection found.")
} else {
  hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  message("Category found.")
}

entrez_col <- if ("ncbi_gene" %in% colnames(hallmark_sets)) "ncbi_gene" else "entrez_gene"

message(sprintf("[ORA] Using '%s' as the Entrez ID column.", entrez_col))

term2gene <- unique(hallmark_sets[, c("gs_name", entrez_col)])
colnames(term2gene) <- c("gs_name", "entrez_gene")
term2gene$entrez_gene <- as.character(term2gene$entrez_gene)

# MSigDB Hallmark enrichment - all DEGs

hallmark_all <- enricher(
  gene = as.character(list_deg),
  TERM2GENE = term2gene,
  universe = as.character(list_universe),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

# Hallmark - upregulated DEGs

hallmark_up <- enricher(
  gene = as.character(list_deg_up),
  TERM2GENE = term2gene,
  universe = as.character(list_universe),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

# Hallmark - downregulated DEGs

hallmark_down <- enricher(
  gene = as.character(list_deg_down),
  TERM2GENE = term2gene,
  universe = as.character(list_universe),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

clean_hallmark_names <- function(x) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(x)
  x@result$Description <- tolower(gsub("_", " ", sub("^HALLMARK_", "", x@result$ID)))
  x
}

hallmark_all  <- clean_hallmark_names(hallmark_all)
hallmark_up   <- clean_hallmark_names(hallmark_up)
hallmark_down <- clean_hallmark_names(hallmark_down)

# Generate enrichment plots

hallmark_all_plot  <- plot_enrichment(hallmark_all,  title = "Hallmark H Enrichment - MDS vs HC (All genes)", color = "#C44E52")
hallmark_up_plot   <- plot_enrichment(hallmark_up,   title = "Hallmark H Enrichment - MDS vs HC (Up-regulated genes)", color = "#009E73")
hallmark_down_plot <- plot_enrichment(hallmark_down, title = "Hallmark H Enrichment - MDS vs HC (Down-regulated genes)", color = "#CC79A7")

# Save enrichment plots

save_plot(hallmark_all_plot, filename = "20-Hallmark_all_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(hallmark_up_plot, filename = "21-Hallmark_up_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)
save_plot(hallmark_down_plot, filename = "22-Hallmark_down_genes.png", dir = paths$plots_gse19429_ora, width = 10, height = 9)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_ORA_GSE19429.txt", dir = paths$logs, label = "Over-Representation Analysis (ORA) - GSE19429")

message("[OUTPUT] Session information saved to:" , paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Over-Representation Analysis completed successfully for GSE19429 (MDS vs Healthy Controls).")
message("=================================================")