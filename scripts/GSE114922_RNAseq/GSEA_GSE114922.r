# =========================================================
# Project      : MDS Genomics - Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE114922 (Illumina HiSeq 4000)
# Samples      : 82 MDS patients + 8 healthy controls, CD34+ bone marrow HSC
# Script       : Gene Set Enrichment Analysis (GSEA)
# Description  : Pre-ranked GSEA using gene-level statistics (e.g. logFC, t-statistics)
#                Functional enrichment against GO, KEGG, Reactome and MSigDB Hallmark
# =========================================================

# =========================================================
#                     Functions
# =========================================================

# Results and lollipop plot 

gsea_lollipop <- function(obj, db_name = "DB", color = "firebrick", title = NULL,
                          top_n = 10, fdr_cutoff = 0.05, min_leading = 10) {
  
  # 1. Standardize input
  
  if ("leadingEdge" %in% names(obj)) {
    
    # fgsea format
    df <- as.data.frame(obj) |>
      dplyr::rename(
        Pathway = pathway,
        p.adjust = padj,
        GeneHit = size
      ) |>
      dplyr::mutate(
        gene = leadingEdge,
        source = "fgsea"
      )
    
  } else {
    
    # clusterProfiler format
    df <- obj@result |>
      dplyr::rename(
        Pathway = Description,
        GeneHit = setSize
      ) |>
      dplyr::mutate(
        gene = strsplit(core_enrichment, "/"),
        source = "clusterProfiler"
      )
  }
  
  # 2. Filter significant
  
  df <- df |>
    dplyr::filter(!is.na(p.adjust), p.adjust < fdr_cutoff)
  
  if (nrow(df) == 0) {
    message("No significant pathways")
    return(NULL)
  }
  
  # 3. Split NES
  
  up <- df |> dplyr::filter(NES > 0)
  down <- df |> dplyr::filter(NES < 0)
  
  # 4. Top selection pathways
  
  select_top <- function(x) {
    x |>
      dplyr::filter(lengths(gene) >= min_leading) |>
      dplyr::arrange(p.adjust, dplyr::desc(GeneHit)) |>
      dplyr::slice_head(n = top_n)
  }
  
  up <- select_top(up)
  down <- select_top(down)
  
  df_top <- dplyr::bind_rows(up, down)
  
  if (nrow(df_top) == 0) {
    message(sprintf("No pathways pass the min_leading filter (>= %d genes)", min_leading))
    return(NULL)
  }
  
  df_expanded <- df_top |>
    tidyr::unnest(gene)
  
  # 6. df for plot - one raw per pathway
  
  df_plot <- df_top |>
    dplyr::mutate(
      FDR = -log10(p.adjust),
      Pathway = reorder(Pathway, NES)
    )
  
  # 7. Lollipop plot
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = NES, y = Pathway)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = NES, yend = Pathway),
      color = "grey70"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = GeneHit, color = FDR)
    ) +
    ggplot2::scale_color_gradient(
      low = "grey60", high = color, name = "-log10(FDR)"
    ) +
    ggplot2::scale_size_continuous(name = "Leading edge size") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = ""
    )
  
  return(list(
    plot = p,
    data = df_expanded,   
    summary = df_top,     
    up = up,
    down = down
  ))
}

# Clean Hallmark names

clean_hallmark_names <- function(x) {
  
  if (is.null(x) || nrow(x) == 0)
    return(x)
  
  x$pathway <- gsub("^HALLMARK_", "", x$pathway)
  x$pathway <- gsub("_", " ", x$pathway)
  x$pathway <- tools::toTitleCase(tolower(x$pathway))
  
  x
}


# =========================================================
#                  Libraries & Setup
# =========================================================

source("SetupEnvironment/00_paths.R")
source("SetupEnvironment/01_environment.R")
source("SetupEnvironment/02_seed.R")
source("SetupEnvironment/03_helper_functions.R")

set_seed(1234)

message("=== Starting GSE114922 GSEA===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Reading differential expression results...")

deg_table <- read.csv(file = file.path(paths$tables, "DEGs_GSE114922.csv"), sep = "," )

message("[LOAD] Reading background gene list...")

background <- read.csv(file = file.path(paths$tables, "background_GSE114922.csv"), sep = "," )

# Ranking 

rank <- data.frame(
  score = deg_table$stat,
  Entrez = deg_table$ENTREZID)

rank_df <- rank |> 
  filter(!is.na(Entrez))

rank_df <- rank_df |> 
  group_by(Entrez) |> 
  slice_max(abs(score), n = 1, with_ties = FALSE) |> 
  ungroup() |> 
  arrange(desc(score))

list_allgenes_ranked <- rank_df$score
names(list_allgenes_ranked) <- rank_df$Entrez

list_allgenes_ranked <- list_allgenes_ranked[is.finite(list_allgenes_ranked)]

summary(list_allgenes_ranked)

# Percentage of genes overlap MDS related genes

mds_like_sets <- msigdbr(species = "Homo sapiens", collection = "C2")

mds_sets <- mds_like_sets |>
  dplyr::filter(grepl("myelodys|myeloid|leuk|hematopoiet", gs_name, ignore.case = TRUE))

mds_genes <- unique(mds_sets$ncbi_gene)

overlap_genes <- intersect(rank_df$Entrez, mds_genes)

percent_overlap <- length(overlap_genes) / length(rank_df$Entrez) * 100

message(sprintf("Overlap MDS gene set: %.2f%% of genes overlap with reference MDS-related signatures", percent_overlap))


# =========================================================
#                  Run GSEA - KEGG Database
# =========================================================

kegg <- gseKEGG(
  geneList = list_allgenes_ranked,
  organism = "hsa",
  keyType = "ncbi-geneid",
  nPerm = 10000,
  pAdjustMethod = "BH",
  verbose = TRUE
)

kegg_df <- as.data.frame(kegg)
res_kegg <- gsea_lollipop(kegg, db_name = "clusterProfiler", color = "#762A83", title = "GSEA enrichment analysis – Reactome pathways")


# =========================================================
#                  Run GSEA - Reactome Database
# =========================================================

reactome <- gsePathway(
  geneList = list_allgenes_ranked,       
  organism = "human",
  pAdjustMethod = "BH",
  minGSSize = 10,               
  maxGSSize = 500,              
  verbose = TRUE
)
res_reactome <- gsea_lollipop(reactome, db_name = "clusterProfiler", color = "darkorange3", title = "GSEA enrichment analysis – Reactome pathways")


# =========================================================
#                  Run GSEA - MSigDB Hallmark 
# =========================================================

message("[GSEA] Retrieving MSigDB Hallmark (H) gene sets...")

hallmark_sets_H <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets_C3 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C3", subcollection = "TFT:GTRD")
hallmark_sets_C7 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C7")
hallmark_sets_C8 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C8")

hm_list_H <- split(
  x = hallmark_sets_H$ncbi_gene, 
  f = hallmark_sets_H$gs_name
)

hm_list_C3 <- split(
  x = hallmark_sets_C3$ncbi_gene, 
  f = hallmark_sets_C3$gs_name
)

hm_list_C7 <- split(
  x = hallmark_sets_C7$ncbi_gene, 
  f = hallmark_sets_C7$gs_name
)

hm_list_C8 <- split(
  x = hallmark_sets_C8$ncbi_gene, 
  f = hallmark_sets_C8$gs_name
)

# Category -> H

hallmarks_H <- fgsea(
  pathways = hm_list_H,
  stats = list_allgenes_ranked
)

## Prepare results for other analyses (leading edge not needed)

hallmarks_H_df <- as.data.frame(hallmarks_H)
hallmarks_H_df_clean <- hallmarks_H_df[, !(names(hallmarks_H_df) %in% "leadingEdge")]

## Process GSEA results and generate lollipop plot

hallmarks_H  <- clean_hallmark_names(hallmarks_H)
res_h <- gsea_lollipop(obj = hallmarks_H, db_name = "Hallmark", color = "firebrick", title = "GSEA enrichment analysis – HallMarks H pathways")

# Category -> C3 (TF targets)

hallmarks_C3 <- fgsea(
  pathways = hm_list_C3,
  stats = list_allgenes_ranked
)
hallmarks_C3 <- clean_hallmark_names(hallmarks_C3)
res_c3 <- gsea_lollipop(obj = hallmarks_C3, db_name = "Hallmark", color = "steelblue", title = "GSEA enrichment analysis – HallMarks C3 pathways")

# Category -> C7 (Immunologic signatures)

hallmarks_C7 <- fgsea(
  pathways = hm_list_C7,
  stats = list_allgenes_ranked
)
hallmarks_C7 <- clean_hallmark_names(hallmarks_C7)
res_c7 <- gsea_lollipop(obj = hallmarks_C7, db_name = "Hallmark", color = "darkgreen", title = "GSEA enrichment analysis – HallMarks C7 pathways")

# Category -> C8 (cell type signatures)

hallmarks_C8 <- fgsea(
  pathways = hm_list_C8,
  stats = list_allgenes_ranked
)
hallmarks_C8 <- clean_hallmark_names(hallmarks_C8)
res_c8 <- gsea_lollipop(obj = hallmarks_C8, db_name = "Hallmark", color = "indianred4", title = "GSEA enrichment analysis – HallMarks C8 pathways")

# -----------------------
# Enrichment score plot
# -----------------------
# IFN-alpha
enrich_IFN_a <- fgsea::plotEnrichment(
  pathway = hm_list_H[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]],
  stats = list_allgenes_ranked) +
  labs(
    title = "Interferon Alpha Response",
    subtitle = "Gene Set Enrichment Analysis",
    x = "Ranked genes",
    y = "Running Enrichment Score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

# IFN-gamma
enrich_IFN_g <- fgsea::plotEnrichment(
  pathway = hm_list_H[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]],
  stats = list_allgenes_ranked) +
  labs(
    title = "Interferon Gamma Response",
    subtitle = "Gene Set Enrichment Analysis",
    x = "Ranked genes",
    y = "Running Enrichment Score"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )


# =========================================================
#             Saving Dataframes  
# =========================================================

# Kegg
save_csv(kegg_df, filename = "results_GSE114922_GSEA_kegg_df.csv", dir = paths$tables_gsea)

# HallMarks H Collection
save_csv(hallmarks_H_df_clean, filename = "results_GSE114922_GSEA_hallmarks_h_df.csv", dir = paths$tables_gsea)


# =========================================================
#             Saving final plots 
# =========================================================

# Reactome 
save_plot(res_reactome[[1]], filename = "23_GSEA_reactome.png", dir = paths$plots_gse114922_gsea, width = 11, height = 9)

# HallMarks 
save_plot(res_h[[1]], filename = "24_GSEA_hallmarks_H.png", dir = paths$plots_gse114922_gsea, width = 11, height = 9)
save_plot(res_c3[[1]], filename = "25_GSEA_hallmarks_C3.png", dir = paths$plots_gse114922_gsea, width = 11, height = 9)
save_plot(res_c7[[1]], filename = "26_GSEA_hallmarks_C7.png", dir = paths$plots_gse114922_gsea, width = 15, height = 9)
save_plot(res_c8[[1]], filename = "27_GSEA_hallmarks_C8.png", dir = paths$plots_gse114922_gsea, width = 11, height = 9)

# Enrichment score
save_plot(enrich_IFN_a, filename = "28_GSEA_enrichmentscore_IFN-A.png", dir = paths$plots_gse114922_gsea, width = 11, height = 9)
save_plot(enrich_IFN_g, filename = "29_GSEA_enrichmentscore_IFN-Gamma.png", dir = paths$plots_gse114922_gsea, width = 11, height = 9)

# Kegg
save_plot(res_kegg[[1]], filename = "30_GSEA_kegg.png", dir = paths$plots_gse114922_gsea, width = 11, height = 9)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_GSEA_GSE114922.txt", dir = paths$logs, label = "Gene Set Enrichment Analysis (GSEA)- GSE114922")

message("[OUTPUT] Session information saved to:" , paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Gene Set Enrichment Analysis completed successfully for GSE114922 (MDS vs Healthy Controls).")
message("=================================================")


